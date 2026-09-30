"""
03_quality_filter.py

Aplica los filtros de calidad a los nanobodies y apunta cuantas filas
sobreviven a cada uno (y cuantas de ellas son de HER2). Los umbrales estan
aqui arriba (y se pueden cambiar desde la terminal, ver abajo) para que
cada decision sea visible.

Filtros, en este orden (mismo orden que la tabla de DECISIONS.md):

  1. Resolucion <= MAX_RES  y  metodo = rayos X o cryo-EM
     (datos de RCSB, calculados en el paso 02). Las estructuras sin datos
     en RCSB se descartan y se cuentan aparte.
  2. Longitud del dominio variable (columna VH) entre 100 y 140 aminoacidos.
     Un dominio VHH tipico mide unos 110-130.

(El filtro "unido a una proteina" y el de falsos nanobodies ya se aplicaron
en el paso 01.)

Como cambiar un umbral sin tocar el codigo (y asi poder comparar):
  python 03_quality_filter.py --max-res 2.5

Salida:
  data/interim/filtered.tsv     -> filas que sobreviven (con 3 columnas extra de RCSB)
  data/interim/filter_log.csv   -> cuantas filas quedan tras cada filtro
"""

import argparse
import csv
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
ENTRADA = BASE / "data" / "raw" / "sabdab_nanobodies.tsv"
PDB_CSV = BASE / "data" / "raw" / "pdb_metadata.csv"
LOG_01 = BASE / "data" / "interim" / "filter_log_01.csv"
SALIDA = BASE / "data" / "interim" / "filtered.tsv"
LOG = BASE / "data" / "interim" / "filter_log.csv"

MAX_RES_POR_DEFECTO = 3.0
METODOS_PERMITIDOS = {"X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"}
LONGITUD_VH = (100, 140)


def limpiar_id(crudo):
    """'pdb_00001n8z' -> '1n8z'"""
    return (crudo or "").strip().split("-")[0][-4:].lower()


def cargar_rcsb(ruta):
    """Devuelve {pdb_id: (metodo, resolucion)} solo para los que tienen datos."""
    datos = {}
    with open(ruta, newline="", encoding="utf-8-sig") as f:
        for fila in csv.DictReader(f):
            try:
                resolucion = float(fila["resolution"]) if fila["resolution"] else None
            except ValueError:
                resolucion = None
            metodo = (fila.get("method") or "").strip().upper()
            if resolucion is not None and metodo:
                datos[fila["pdb_id"]] = (metodo, resolucion)
    return datos


def leer_log_01():
    """Recupera las lineas del log del paso 01, si existe."""
    if not LOG_01.exists():
        return []
    with open(LOG_01, newline="", encoding="utf-8") as f:
        lector = csv.reader(f)
        next(lector, None)
        return [(fila[0], int(fila[1]), int(fila[2])) for fila in lector if len(fila) == 3]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-res", type=float, default=MAX_RES_POR_DEFECTO,
                        help="resolucion maxima en A (por defecto 3.0)")
    args = parser.parse_args()
    max_res = args.max_res

    for ruta in (ENTRADA, PDB_CSV):
        if not ruta.exists():
            sys.exit(f"No encuentro {ruta}. Corre antes los scripts 01 y 02.")

    rcsb = cargar_rcsb(PDB_CSV)

    with open(ENTRADA, newline="", encoding="utf-8-sig", errors="replace") as f:
        lector = csv.DictReader(f, delimiter="\t")
        columnas = lector.fieldnames or []
        filas = list(lector)

    por_minuscula = {c.lower(): c for c in columnas}
    col_pdb = por_minuscula.get("pdb") or por_minuscula.get("instance")
    col_vh = por_minuscula.get("vh")
    col_her2 = por_minuscula.get("is_her2")
    for nombre, col in [("PDB", col_pdb), ("VH", col_vh), ("is_her2", col_her2)]:
        if col is None:
            sys.exit(f"Falta la columna {nombre}. Corre otra vez el script 01. Columnas: {columnas}")

    def n_her2(lista):
        return sum(1 for f in lista if f[col_her2] == "True")

    log = leer_log_01()

    # Filtro 1: resolucion + metodo (datos de RCSB)
    paso1 = []
    sin_datos = 0
    for fila in filas:
        pdb_id = limpiar_id(fila[col_pdb])
        if pdb_id not in rcsb:
            sin_datos += 1
            continue
        metodo, resolucion = rcsb[pdb_id]
        if resolucion <= max_res and metodo in METODOS_PERMITIDOS:
            fila = dict(fila)
            fila["pdb_id_4"] = pdb_id
            fila["rcsb_method"] = metodo
            fila["rcsb_resolution"] = resolucion
            paso1.append(fila)
    log.append((f"resolucion_y_metodo (<={max_res}A, rayosX/cryoEM)", len(paso1), n_her2(paso1)))

    # Filtro 2: longitud del dominio variable
    paso2 = [f for f in paso1
             if LONGITUD_VH[0] <= len((f[col_vh] or "").strip()) <= LONGITUD_VH[1]]
    log.append((f"longitud_dominio_vhh ({LONGITUD_VH[0]}-{LONGITUD_VH[1]}aa)", len(paso2), n_her2(paso2)))

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    columnas_salida = columnas + ["pdb_id_4", "rcsb_method", "rcsb_resolution"]
    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        escritor = csv.DictWriter(f, fieldnames=columnas_salida, delimiter="\t")
        escritor.writeheader()
        escritor.writerows(paso2)

    with open(LOG, "w", newline="", encoding="utf-8") as f:
        escritor = csv.writer(f)
        escritor.writerow(["paso", "n_entradas", "n_her2"])
        escritor.writerows(log)

    print(f"{'Paso':<75}{'filas':>8}{'de HER2':>9}")
    for nombre, n, h in log:
        print(f"{nombre:<75}{n:>8}{h:>9}")
    if sin_datos:
        print(f"\nOjo: {sin_datos} filas se descartaron porque RCSB no dio resolucion o metodo "
              "para su estructura (por ejemplo, estructuras de NMR o IDs que RCSB no reconoce).")
    print(f"\nGuardado: {SALIDA}")
    print(f"Guardado: {LOG}")


if __name__ == "__main__":
    main()
