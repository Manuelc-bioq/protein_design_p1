"""
05_build_release.py

Empaqueta el dataset final en una carpeta con version:

  data/releases/<version>/
      sequences.fasta          secuencias de los nanobodies representantes
      dataset.tsv              tabla completa de esos representantes
      cluster_assignments.tsv  a que grupo pertenecia cada fila antes de deduplicar
      filter_log.csv           cuantas filas sobrevivieron a cada paso (con los umbrales)
      MANIFEST.txt             fecha, numero de secuencias, huellas SHA-256 y tu nota

Nunca sobrescribe una version anterior: si la carpeta ya existe, se niega.

Uso:
  python 05_build_release.py --version v1.0 --note "Primera version, HER2 nanobodies"
"""

import argparse
import hashlib
import shutil
import sys
from datetime import date
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
FASTA = BASE / "data" / "processed" / "dedup_representatives.fasta"
TABLA = BASE / "data" / "processed" / "dedup_representatives.tsv"
CLUSTERS = BASE / "data" / "processed" / "cluster_assignments.tsv"
LOG = BASE / "data" / "interim" / "filter_log.csv"
CHANGELOG = BASE / "CHANGELOG.md"


def sha256(ruta):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(8192), b""):
            h.update(bloque)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", required=True)
    parser.add_argument("--note", required=True)
    args = parser.parse_args()

    for ruta in (FASTA, TABLA, CLUSTERS, LOG):
        if not ruta.exists():
            sys.exit(f"Falta {ruta}. Corre antes los scripts 01 a 04.")

    carpeta = BASE / "data" / "releases" / args.version
    if carpeta.exists():
        sys.exit(f"La version {args.version} ya existe ({carpeta}). "
                 "Elige otro numero (por ejemplo v1.1): las versiones no se sobrescriben.")
    carpeta.mkdir(parents=True)

    shutil.copy(FASTA, carpeta / "sequences.fasta")
    shutil.copy(TABLA, carpeta / "dataset.tsv")
    shutil.copy(CLUSTERS, carpeta / "cluster_assignments.tsv")
    shutil.copy(LOG, carpeta / "filter_log.csv")

    with open(carpeta / "sequences.fasta", encoding="utf-8") as f:
        n_secuencias = sum(1 for linea in f if linea.startswith(">"))
    huella_fasta = sha256(carpeta / "sequences.fasta")
    huella_tabla = sha256(carpeta / "dataset.tsv")

    with open(carpeta / "MANIFEST.txt", "w", encoding="utf-8") as f:
        f.write(f"version: {args.version}\n")
        f.write(f"fecha: {date.today().isoformat()}\n")
        f.write("objetivo: todos los nanobodies (VHH) unidos a una proteina; HER2 marcado (is_her2)\n")
        f.write(f"n_secuencias: {n_secuencias}\n")
        f.write(f"sha256(sequences.fasta): {huella_fasta}\n")
        f.write(f"sha256(dataset.tsv): {huella_tabla}\n")
        f.write(f"nota: {args.note}\n")

    with open(CHANGELOG, "a", encoding="utf-8") as f:
        f.write(f"\n## {args.version} - {date.today().isoformat()}\n"
                f"- {args.note}\n"
                f"- n_secuencias: {n_secuencias}\n"
                f"- sha256(sequences.fasta): {huella_fasta}\n")

    print(f"Version {args.version} creada en {carpeta}")
    print(f"  {n_secuencias} secuencias | sha256 = {huella_fasta[:12]}...")


if __name__ == "__main__":
    main()
