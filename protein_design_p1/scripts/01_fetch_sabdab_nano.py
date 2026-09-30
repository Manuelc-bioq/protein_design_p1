"""
01_fetch_sabdab_nano.py

NO DESCARGA NADA. Lee la tabla completa de SAbDab que bajaste a mano
(data/raw/sabdab_nano_all.tsv) y se queda con TODOS los nanobodies que estan
unidos a una proteina. Los que son de HER2 se marcan con una columna nueva
(is_her2 = True), pero NO se separan del resto.

Por que todos y no solo HER2:
  En el PDB hay muy pocos nanobodies anti-HER2 validos (3 secuencias
  utilizables). Con tan pocos, deduplicar y hacer estadistica no tiene
  sentido. Con todos los nanobodies (~900 secuencias distintas) si, y HER2
  sigue identificado para P3 y P4.

Pasos, en este orden (cada uno apunta cuantas filas quedan):
  1. tabla_completa_sabdab             todas las filas de la tabla
  2. nanobody (type = SD-H)            una sola cadena pesada, sin cadena ligera
  3. antigeno_proteina                 el antigeno es una proteina y hay cadena de antigeno
  4. sin_falsos_nanobodies             se quitan las filas cuya CDR-H3 es identica a la de
                                       un anticuerpo convencional de la misma tabla
                                       (son cadenas pesadas de anticuerpos normales
                                       mal etiquetadas, como 8ffj = trastuzumab)

Nota: SAbDab etiqueta los nanobodies como "SD-H" (single domain, heavy).
"SD-L" y "VNAR" (tiburon) NO son nanobodies de camelido y no entran.

Salida:
  data/raw/sabdab_nanobodies.tsv      -> filas que pasan (con la columna is_her2)
  data/interim/filter_log_01.csv      -> cuantas filas quedan en cada paso (y cuantas de HER2)
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
ENTRADA = BASE / "data" / "raw" / "sabdab_nano_all.tsv"
SALIDA = BASE / "data" / "raw" / "sabdab_nanobodies.tsv"
LOG = BASE / "data" / "interim" / "filter_log_01.csv"

# Formas en las que puede aparecer HER2 en el nombre del antigeno
HER2_KEYWORDS = ["her2", "her-2", "erbb2", "erbb-2", "neu proto-oncogene"]

# Valores que significan "aqui no hay nada"
VACIOS = {"", "NA", "N/A", "NAN", "NONE", "NULL", "-"}

# Una CDR-H3 mas corta que esto podria coincidir por casualidad; no se usa para comparar
MIN_LONGITUD_CDR3 = 8


def vacio(valor):
    return (valor or "").strip().upper() in VACIOS


def leer_tabla(ruta):
    """Lee un CSV o TSV y detecta solo el separador. Devuelve (columnas, filas)."""
    with open(ruta, newline="", encoding="utf-8-sig", errors="replace") as f:
        primera_linea = f.readline()
        f.seek(0)
        candidatos = {"\t": primera_linea.count("\t"),
                      ",": primera_linea.count(","),
                      ";": primera_linea.count(";")}
        separador = max(candidatos, key=candidatos.get)
        lector = csv.DictReader(f, delimiter=separador)
        filas = list(lector)
        return lector.fieldnames or [], filas


def buscar_columna(columnas, *nombres):
    """Busca una columna sin importar mayusculas/minusculas."""
    por_minuscula = {c.lower(): c for c in columnas}
    for nombre in nombres:
        if nombre.lower() in por_minuscula:
            return por_minuscula[nombre.lower()]
    return None


def es_her2(fila, col_antigeno):
    nombre = (fila.get(col_antigeno) or "").lower()
    return any(k in nombre for k in HER2_KEYWORDS)


def main():
    parser = argparse.ArgumentParser()
    # Se acepta (y se ignora) por si escribes el comando antiguo con este flag
    parser.add_argument("--skip-download", action="store_true", help=argparse.SUPPRESS)
    parser.parse_args()

    if not ENTRADA.exists():
        sys.exit(f"No encuentro {ENTRADA}\n"
                 "Copia ahi la tabla completa que bajaste de la web de SAbDab.")

    with open(ENTRADA, "rb") as f:
        inicio = f.read(200).lstrip().lower()
    if inicio.startswith(b"<!doctype") or inicio.startswith(b"<html"):
        sys.exit("El archivo sabdab_nano_all.tsv es una pagina web (HTML), no la tabla.\n"
                 "Vuelve a copiar aqui el archivo grande (unos 17 MB) que bajaste de la web.")

    columnas, filas = leer_tabla(ENTRADA)

    col_antigeno = buscar_columna(columnas, "antigen_name")
    col_ag_tipo = buscar_columna(columnas, "antigen_type")
    col_ag_cadena = buscar_columna(columnas, "antigen_chain")
    col_lchain = buscar_columna(columnas, "Lchain")
    col_vl = buscar_columna(columnas, "VL")
    col_vh = buscar_columna(columnas, "VH")
    col_cdr3 = buscar_columna(columnas, "CDR-H3")
    col_tipo = buscar_columna(columnas, "type")

    obligatorias = [("antigen_name", col_antigeno), ("antigen_type", col_ag_tipo),
                    ("antigen_chain", col_ag_cadena), ("VH", col_vh), ("type", col_tipo)]
    faltan = [n for n, c in obligatorias if c is None]
    if faltan:
        sys.exit(f"Faltan columnas necesarias: {faltan}\n"
                 f"Columnas que hay en tu tabla: {columnas}")

    # Paso 2: nanobody = type SD-H
    nanos = [f for f in filas if (f.get(col_tipo) or "").strip().upper() == "SD-H"]

    # Paso 3: unido a una proteina
    con_proteina = [f for f in nanos
                    if (f.get(col_ag_tipo) or "").strip().upper() == "PROTEIN"
                    and not vacio(f.get(col_ag_cadena))]

    # Paso 4: fuera los "falsos nanobodies"
    #   CDR-H3 de los anticuerpos convencionales (los que SI tienen cadena ligera)
    cdr3_convencionales = set()
    if col_cdr3 and col_vl:
        for f in filas:
            cdr3 = (f.get(col_cdr3) or "").strip()
            if (f.get(col_tipo) or "").strip().upper() != "SD-H" and not vacio(f.get(col_vl)) \
                    and len(cdr3) >= MIN_LONGITUD_CDR3:
                cdr3_convencionales.add(cdr3)
    validos, falsos = [], []
    for f in con_proteina:
        cdr3 = (f.get(col_cdr3) or "").strip() if col_cdr3 else ""
        (falsos if cdr3 in cdr3_convencionales else validos).append(f)

    for f in validos:
        f["is_her2"] = str(es_her2(f, col_antigeno))

    def n_her2(lista):
        return sum(1 for f in lista if es_her2(f, col_antigeno))

    log = [
        ("tabla_completa_sabdab", len(filas), n_her2(filas)),
        ("nanobody_SD-H (sin cadena ligera)", len(nanos), n_her2(nanos)),
        ("antigeno_proteina", len(con_proteina), n_her2(con_proteina)),
        ("sin_falsos_nanobodies (CDR-H3 distinta de anticuerpos convencionales)",
         len(validos), n_her2(validos)),
    ]

    print(f"{'Paso':<75}{'filas':>8}{'de HER2':>9}")
    for nombre, n, h in log:
        print(f"{nombre:<75}{n:>8}{h:>9}")

    print("\nTipos de molecula en la tabla completa (columna 'type'):")
    for tipo, n in Counter((f.get(col_tipo) or "").strip() for f in filas).most_common():
        print(f"  {n:6d} x {tipo}")

    print(f"\nFilas quitadas por parecer falsos nanobodies: {len(falsos)}")
    for f in falsos:
        if es_her2(f, col_antigeno):
            print(f"  (de HER2) {f.get('PDB', '')}  CDR-H3 = {f.get(col_cdr3)}")

    print(f"\nEstructuras PDB distintas en el resultado: "
          f"{len({(f.get('PDB') or '').split('-')[0] for f in validos})}")
    print(f"Secuencias VH distintas (100% identicas): {len({f[col_vh] for f in validos})}")

    if not validos:
        sys.exit("\nNo ha quedado ninguna fila. Pasame la tabla de arriba.")

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        escritor = csv.DictWriter(f, fieldnames=columnas + ["is_her2"], delimiter="\t")
        escritor.writeheader()
        escritor.writerows(validos)

    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "w", newline="", encoding="utf-8") as f:
        escritor = csv.writer(f)
        escritor.writerow(["paso", "n_entradas", "n_her2"])
        escritor.writerows(log)

    print(f"\nGuardado: {SALIDA}")
    print(f"Guardado: {LOG}")


if __name__ == "__main__":
    main()
