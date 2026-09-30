"""
02_fetch_pdb_metadata.py

Para cada estructura del PDB que aparece en data/raw/sabdab_nanobodies.tsv,
pregunta a RCSB (la fuente oficial del PDB) por el metodo experimental y la
resolucion. Necesita internet.

Son unas 1.800 estructuras: tarda unos minutos. Hace 4 consultas a la vez.

Se puede parar (Ctrl+C) y volver a lanzar: lo que ya se descargo bien se
guarda y no se repite; solo se reintentan las que fallaron.

Por que RCSB y no la columna "resolution" de SAbDab:
  SAbDab es una copia con anotaciones; RCSB es la fuente original. Si hay
  discrepancia, manda RCSB.

Nota sobre los identificadores:
  SAbDab escribe el codigo PDB como "pdb_00001n8z". El codigo real de PDB
  son solo los 4 ultimos caracteres de esa parte: "1n8z".

Salida:
  data/raw/pdb_metadata.csv   (una fila por estructura PDB unica)
"""

import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

BASE = Path(__file__).resolve().parents[1]
ENTRADA = BASE / "data" / "raw" / "sabdab_nanobodies.tsv"
SALIDA = BASE / "data" / "raw" / "pdb_metadata.csv"
URL_RCSB = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
COLUMNAS = ["pdb_id", "method", "resolution", "deposit_date", "error"]
CONSULTAS_A_LA_VEZ = 4


def limpiar_id(crudo):
    """'pdb_00001n8z' (o 'pdb_00001n8z-B-A') -> '1n8z'"""
    primera_parte = (crudo or "").strip().split("-")[0]
    return primera_parte[-4:].lower()


def leer_ids(ruta):
    with open(ruta, newline="", encoding="utf-8-sig", errors="replace") as f:
        lector = csv.DictReader(f, delimiter="\t")
        por_minuscula = {c.lower(): c for c in (lector.fieldnames or [])}
        col = por_minuscula.get("pdb") or por_minuscula.get("instance")
        if col is None:
            sys.exit(f"No encuentro la columna PDB. Columnas: {lector.fieldnames}")
        return sorted({limpiar_id(fila[col]) for fila in lector if fila.get(col)})


def leer_ya_descargados(ruta):
    """Devuelve {pdb_id: fila} con las estructuras que ya se consultaron SIN error."""
    hechos = {}
    if ruta.exists():
        with open(ruta, newline="", encoding="utf-8") as f:
            for fila in csv.DictReader(f):
                if not fila.get("error"):
                    hechos[fila["pdb_id"]] = fila
    return hechos


def pedir_a_rcsb(pdb_id):
    vacio = {"pdb_id": pdb_id, "method": None, "resolution": None, "deposit_date": None}
    respuesta = None
    error = None
    for _ in range(3):  # dos reintentos si falla la conexion
        try:
            respuesta = requests.get(URL_RCSB.format(pdb_id=pdb_id), timeout=30)
            break
        except requests.RequestException as e:
            error = type(e).__name__
            time.sleep(1)
    if respuesta is None:
        return {**vacio, "error": error}
    if respuesta.status_code != 200:
        return {**vacio, "error": f"HTTP {respuesta.status_code}"}

    datos = respuesta.json()
    try:
        metodo = datos["exptl"][0]["method"]
    except (KeyError, IndexError, TypeError):
        metodo = None
    try:
        resolucion = datos["rcsb_entry_info"]["resolution_combined"][0]
    except (KeyError, IndexError, TypeError):
        resolucion = None
    return {
        "pdb_id": pdb_id,
        "method": metodo,
        "resolution": resolucion,
        "deposit_date": (datos.get("rcsb_accession_info") or {}).get("deposit_date"),
        "error": None,
    }


def guardar(ruta, filas_por_id):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        escritor = csv.DictWriter(f, fieldnames=COLUMNAS)
        escritor.writeheader()
        for pdb_id in sorted(filas_por_id):
            escritor.writerow({c: filas_por_id[pdb_id].get(c) for c in COLUMNAS})


def main():
    if not ENTRADA.exists():
        sys.exit(f"No encuentro {ENTRADA}. Corre antes el script 01.")
    ids = leer_ids(ENTRADA)
    filas = leer_ya_descargados(SALIDA)
    filas = {i: f for i, f in filas.items() if i in set(ids)}
    pendientes = [i for i in ids if i not in filas]
    print(f"{len(ids)} estructuras PDB unicas | ya descargadas: {len(filas)} | "
          f"por consultar: {len(pendientes)}")

    try:
        with ThreadPoolExecutor(max_workers=CONSULTAS_A_LA_VEZ) as pool:
            trabajos = [pool.submit(pedir_a_rcsb, i) for i in pendientes]
            for n, trabajo in enumerate(as_completed(trabajos), 1):
                resultado = trabajo.result()
                filas[resultado["pdb_id"]] = resultado
                if n % 100 == 0:
                    print(f"  ... {n}/{len(pendientes)}")
                    guardar(SALIDA, filas)  # punto de guardado
    except KeyboardInterrupt:
        print("\nInterrumpido. Guardo lo que hay; vuelve a lanzar el script para continuar.")
        guardar(SALIDA, filas)
        sys.exit(1)

    guardar(SALIDA, filas)
    con_datos = sum(1 for f in filas.values() if f["resolution"] not in (None, ""))
    con_error = sum(1 for f in filas.values() if f["error"])
    print(f"Con resolucion disponible: {con_datos} | con error: {con_error}")
    if con_error:
        print("Hay estructuras con error: vuelve a lanzar el script y se reintentan solo esas.")
    print(f"Guardado: {SALIDA}")


if __name__ == "__main__":
    main()
