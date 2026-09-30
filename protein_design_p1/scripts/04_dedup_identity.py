"""
04_dedup_identity.py

Deduplica por identidad de secuencia SIN instalar nada (Python puro).
Sustituye a CD-HIT, que es dificil de instalar en Windows.

Como funciona, paso a paso:
  1. Agrupa las filas con la secuencia VH EXACTAMENTE igual (muchas estructuras
     del PDB repiten el mismo nanobody). Queda una lista de secuencias unicas.
  2. Ordena esas secuencias de mejor a peor: menor resolucion primero (mas
     nitida); a igualdad, la mas larga.
  3. Recorre la lista. Cada secuencia se compara con los "representantes" que
     ya existen:
       - si se parece a alguno en >= IDENTIDAD (por defecto 95%), entra en el
         grupo del que MAS se parece;
       - si no se parece a ninguno, crea un grupo nuevo y ella es su
         representante.
  4. Como la lista va de mejor a peor, el representante de cada grupo es
     siempre la estructura de mejor resolucion.

Que es "identidad":
  Se alinean las dos secuencias (Needleman-Wunsch: se colocan una debajo de
  otra, permitiendo huecos) y se cuentan las posiciones con el mismo
  aminoacido, divididas entre la longitud de la secuencia MAS CORTA
  (la misma definicion que usa CD-HIT).

Para no tardar media hora, antes de alinear se descartan los pares que
comparten muy pocos "trozos de 4 aminoacidos" (k-mers): dos secuencias con
>= 95% de identidad comparten casi todos, asi que esos pares no pueden
llegar al umbral. Es un atajo; el resultado es el mismo que comparar todo
con todo (se comprobo con los datos reales, ver DECISIONS.md).

Como cambiar el umbral y comparar:
  python 04_dedup_identity.py --identity 0.90

Ojo (limitacion conocida, ver DECISIONS.md): dos nanobodies con el mismo
armazon pero CDR-H3 distinta pueden acabar en el mismo grupo. Por eso las
tablas de salida guardan la CDR-H3 y este script cuenta los grupos donde no
coinciden.

Salida:
  data/processed/dedup_representatives.fasta
  data/processed/dedup_representatives.tsv    (filas completas de los representantes)
  data/processed/cluster_assignments.tsv      (a que grupo pertenece cada fila)
  y una linea nueva en data/interim/filter_log.csv
"""

import argparse
import csv
import math
import sys
import time
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
ENTRADA = BASE / "data" / "interim" / "filtered.tsv"
LOG = BASE / "data" / "interim" / "filter_log.csv"
SAL_FASTA = BASE / "data" / "processed" / "dedup_representatives.fasta"
SAL_TSV = BASE / "data" / "processed" / "dedup_representatives.tsv"
SAL_CLUSTERS = BASE / "data" / "processed" / "cluster_assignments.tsv"

IDENTIDAD_POR_DEFECTO = 0.95
K = 4  # tamano de los "trozos" del prefiltro


def identidad(a, b):
    """Identidad de secuencia = coincidencias / longitud de la mas corta."""
    if a == b:
        return 1.0
    n, m = len(a), len(b)
    ACIERTO, FALLO, HUECO = 2, -1, -2
    puntos = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        puntos[i][0] = i * HUECO
    for j in range(1, m + 1):
        puntos[0][j] = j * HUECO
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = puntos[i - 1][j - 1] + (ACIERTO if a[i - 1] == b[j - 1] else FALLO)
            puntos[i][j] = max(diagonal, puntos[i - 1][j] + HUECO, puntos[i][j - 1] + HUECO)
    # Recorrido hacia atras para contar coincidencias
    i, j, coincidencias = n, m, 0
    while i > 0 and j > 0:
        igual = a[i - 1] == b[j - 1]
        if puntos[i][j] == puntos[i - 1][j - 1] + (ACIERTO if igual else FALLO):
            coincidencias += 1 if igual else 0
            i -= 1
            j -= 1
        elif puntos[i][j] == puntos[i - 1][j] + HUECO:
            i -= 1
        else:
            j -= 1
    return coincidencias / min(n, m)


def trozos(secuencia):
    """Todos los trozos de K aminoacidos consecutivos (en orden, con repeticiones)."""
    return [secuencia[i:i + K] for i in range(len(secuencia) - K + 1)]


def puede_llegar_al_umbral(trozos_corta, conjunto_larga, largo_corta, umbral):
    """Cota barata: si comparten muy pocos trozos, la identidad no llega al umbral.

    Cada aminoacido que NO coincide destruye como mucho K trozos de la secuencia corta.
    Con identidad >= umbral hay como mucho (1 - umbral) * largo aminoacidos distintos.
    Se deja un margen extra por si hay huecos (indels) en el alineamiento.
    """
    distintos_max = largo_corta - math.ceil(umbral * largo_corta)
    minimo = len(trozos_corta) - K * distintos_max - 4 * K
    if minimo <= 0:
        return True
    compartidos = sum(1 for t in trozos_corta if t in conjunto_larga)
    return compartidos >= minimo


def agrupar(secuencias, umbral, usar_prefiltro=True):
    """secuencias: lista ya ordenada de mejor a peor (unicas).
    Devuelve (grupo_de[i], identidad_al_representante[i], indices_representantes)."""
    trozos_de = [trozos(s) for s in secuencias]
    conjuntos = [set(t) for t in trozos_de]
    representantes = []          # indices de secuencias que son representantes
    grupo_de, ident_de = [], []
    for i, secuencia in enumerate(secuencias):
        mejor_grupo, mejor_id = None, 0.0
        for id_grupo, r in enumerate(representantes):
            rep = secuencias[r]
            if len(secuencia) <= len(rep):
                corta_trozos, larga_conjunto, largo = trozos_de[i], conjuntos[r], len(secuencia)
            else:
                corta_trozos, larga_conjunto, largo = trozos_de[r], conjuntos[i], len(rep)
            if usar_prefiltro and not puede_llegar_al_umbral(corta_trozos, larga_conjunto, largo, umbral):
                continue
            valor = identidad(secuencia, rep)
            if valor > mejor_id:
                mejor_grupo, mejor_id = id_grupo, valor
        if mejor_grupo is not None and mejor_id >= umbral:
            grupo_de.append(mejor_grupo)
            ident_de.append(mejor_id)
        else:
            representantes.append(i)
            grupo_de.append(len(representantes) - 1)
            ident_de.append(1.0)
    return grupo_de, ident_de, representantes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--identity", type=float, default=IDENTIDAD_POR_DEFECTO,
                        help="identidad minima para considerar dos secuencias 'la misma' (0-1)")
    parser.add_argument("--sin-prefiltro", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    umbral = args.identity

    if not ENTRADA.exists():
        sys.exit(f"No encuentro {ENTRADA}. Corre antes el script 03.")
    with open(ENTRADA, newline="", encoding="utf-8-sig") as f:
        lector = csv.DictReader(f, delimiter="\t")
        columnas = lector.fieldnames or []
        filas = list(lector)
    if not filas:
        sys.exit("El filtro de calidad no dejo ninguna fila. Mira el log del paso 03.")

    por_minuscula = {c.lower(): c for c in columnas}
    col_id = por_minuscula["instance"]
    col_vh = por_minuscula["vh"]
    col_cdr3 = por_minuscula.get("cdr-h3")
    col_her2 = por_minuscula["is_her2"]

    # Paso 1: secuencias exactamente iguales -> una sola secuencia unica
    filas.sort(key=lambda f: (float(f["rcsb_resolution"]), -len(f[col_vh]), f[col_id]))
    filas_de_secuencia = {}      # secuencia -> filas (la primera es la mejor)
    for fila in filas:
        filas_de_secuencia.setdefault(fila[col_vh].strip(), []).append(fila)
    secuencias = list(filas_de_secuencia)   # ya en orden de mejor a peor
    print(f"{len(filas)} filas -> {len(secuencias)} secuencias VH distintas (100% identicas agrupadas)")

    # Pasos 2-4: agrupar secuencias unicas por identidad
    inicio = time.time()
    grupo_de, ident_de, reps = agrupar(secuencias, umbral, usar_prefiltro=not args.sin_prefiltro)
    print(f"Agrupado en {time.time() - inicio:.1f} s")

    # Reparto de filas por grupo
    filas_de_grupo = {}
    asignaciones = []            # (fila, grupo, identidad, es_representante)
    for indice, secuencia in enumerate(secuencias):
        for j, fila in enumerate(filas_de_secuencia[secuencia]):
            g = grupo_de[indice]
            es_rep = (indice == reps[g]) and j == 0
            asignaciones.append((fila, g, ident_de[indice], es_rep))
            filas_de_grupo.setdefault(g, []).append(fila)

    representante_de = {g: filas_de_secuencia[secuencias[r]][0] for g, r in enumerate(reps)}

    def resumen_grupo(g):
        miembros = filas_de_grupo[g]
        return {
            "cluster_size": len(miembros),
            "cluster_n_sequences": len({m[col_vh].strip() for m in miembros}),
            "cluster_has_her2": any(m[col_her2] == "True" for m in miembros),
            "cluster_n_cdr_h3": len({m[col_cdr3] for m in miembros}) if col_cdr3 else "",
        }

    SAL_FASTA.parent.mkdir(parents=True, exist_ok=True)
    with open(SAL_FASTA, "w", encoding="utf-8") as f:
        for g in range(len(reps)):
            rep = representante_de[g]
            f.write(f">{rep[col_id]} pdb={rep['pdb_id_4']} res={rep['rcsb_resolution']} "
                    f"grupo={g} her2={rep[col_her2]}\n{rep[col_vh].strip()}\n")

    extra = ["cluster_id", "cluster_size", "cluster_n_sequences", "cluster_has_her2", "cluster_n_cdr_h3"]
    with open(SAL_TSV, "w", newline="", encoding="utf-8") as f:
        escritor = csv.DictWriter(f, fieldnames=columnas + extra, delimiter="\t")
        escritor.writeheader()
        for g in range(len(reps)):
            escritor.writerow({**representante_de[g], "cluster_id": g, **resumen_grupo(g)})

    with open(SAL_CLUSTERS, "w", newline="", encoding="utf-8") as f:
        escritor = csv.writer(f, delimiter="\t")
        escritor.writerow(["entry", "pdb_id", "cluster_id", "is_representative",
                           "identity_to_representative", "cdr_h3", "resolution", "is_her2"])
        for fila, g, ident, es_rep in asignaciones:
            escritor.writerow([fila[col_id], fila["pdb_id_4"], g, es_rep, f"{ident:.3f}",
                               fila[col_cdr3] if col_cdr3 else "", fila["rcsb_resolution"],
                               fila[col_her2]])

    # Anadir una linea al log (sustituyendo la de una ejecucion anterior)
    n_her2_grupos = sum(1 for g in range(len(reps)) if resumen_grupo(g)["cluster_has_her2"])
    nombre_paso = f"tras_deduplicacion ({umbral * 100:.0f}% identidad)"
    lineas = []
    if LOG.exists():
        with open(LOG, newline="", encoding="utf-8") as f:
            lineas = [fila for fila in csv.reader(f) if fila and not fila[0].startswith("tras_deduplicacion")]
    else:
        lineas = [["paso", "n_entradas", "n_her2"]]
    lineas.append([nombre_paso, str(len(reps)), str(n_her2_grupos)])
    with open(LOG, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(lineas)

    # Resumen por pantalla
    print(f"\nUmbral de identidad: {umbral * 100:.0f}%")
    print(f"Filas antes de deduplicar:      {len(filas)}")
    print(f"Grupos (secuencias distintas):  {len(reps)}")
    print(f"Grupos que incluyen HER2:       {n_her2_grupos}")
    print(f"Representantes que son HER2:    {sum(1 for g in range(len(reps)) if representante_de[g][col_her2] == 'True')}")

    grandes = sorted(range(len(reps)), key=lambda g: -len(filas_de_grupo[g]))[:5]
    print("\nLos 5 grupos mas grandes (filas | secuencias distintas | CDR-H3 distintas):")
    for g in grandes:
        r = resumen_grupo(g)
        print(f"  grupo {g}: {r['cluster_size']} filas | {r['cluster_n_sequences']} secuencias | "
              f"{r['cluster_n_cdr_h3']} CDR-H3 | representante {representante_de[g]['pdb_id_4']}")

    if col_cdr3:
        mezclados = [g for g in range(len(reps)) if resumen_grupo(g)["cluster_n_cdr_h3"] > 1]
        print(f"\nGrupos con MAS DE UNA CDR-H3 distinta: {len(mezclados)} de {len(reps)} "
              "(limitacion conocida, ver DECISIONS.md apartado 3)")
        for g in mezclados[:5]:
            cdr3 = sorted({m[col_cdr3] for m in filas_de_grupo[g]})
            print(f"  grupo {g}: {', '.join(cdr3[:4])}{' ...' if len(cdr3) > 4 else ''}")

    print(f"\nGuardado: {SAL_FASTA}")
    print(f"Guardado: {SAL_TSV}")
    print(f"Guardado: {SAL_CLUSTERS}")


if __name__ == "__main__":
    main()
