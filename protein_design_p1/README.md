# Nanobody–protein complex dataset (PDB, curated and deduplicated)

A versioned dataset of **551 non-redundant nanobody (VHH) sequences** taken from
experimentally solved nanobody–protein complexes in the PDB, filtered by structural
quality and deduplicated at 95 % sequence identity. Nanobodies against **HER2/ERBB2**
are flagged (`is_her2`) so they can be pulled out as a subset.

This is project **P1** of a self-directed learning roadmap in computational protein
design. It is the data foundation for later projects (train/test splitting, interface
analysis, comparing designed binders against natural ones).

## What is in the release

`data/releases/v1.0/`

| File | Content |
|---|---|
| `sequences.fasta` | 551 representative VHH sequences (one per 95 % identity cluster) |
| `dataset.tsv` | Full SAbDab annotation for each representative, plus RCSB method/resolution, `is_her2`, and cluster summary columns |
| `cluster_assignments.tsv` | Every structure that passed the filters (1,791 rows) and the cluster it was assigned to |
| `filter_log.csv` | How many rows survived each step, and how many of them are HER2 |
| `MANIFEST.txt` | Date, count, SHA-256 checksums of the FASTA and TSV |

## How the data was filtered

| Step | Rows | of which HER2 |
|---|---:|---:|
| Full SAbDab table | 22,264 | 32 |
| Nanobody (`type = SD-H`, no light chain) | 4,614 | 5 |
| Bound to a protein antigen | 3,212 | 5 |
| Removed mislabelled heavy chains (CDR-H3 identical to a conventional antibody in the same table, e.g. 8ffj = trastuzumab heavy chain) | 3,202 | 4 |
| Resolution ≤ 3.0 Å, X-ray or cryo-EM (values from RCSB) | 1,797 | 3 |
| VH domain length 100–140 aa | 1,791 | 3 |
| **Deduplicated at 95 % identity (best-resolution representative)** | **551** | **3** |

The reasoning behind every threshold is in [`DECISIONS.md`](DECISIONS.md) (in Spanish).

## How this was built — and what I did

The code was written with **Claude (Anthropic) as an AI coding assistant**. I am a
second-year biochemistry student and do not yet write Python at this level. My part was:

- choosing the target and deciding the scope, using the data: the original plan was
  anti-HER2 nanobodies only, but the PDB contains just 3 usable sequences, too few to
  deduplicate or describe, so I widened it to all nanobody–protein complexes with HER2 flagged;
- deciding and documenting the curation thresholds (resolution, method, length, identity);
- running every step and checking the output. For example, I questioned a drop from
  32 to 24 structures between two steps, which exposed how SAbDab IDs map to PDB entries,
  and later the finding that the first "HER2 nanobodies" were in fact conventional antibodies.

What I did **not** do: write the scripts from scratch, or have the results checked by an
independent reviewer.

## Known limitations

- **HER2 subset is tiny**: 3 sequences from 2 PDB entries (5my6, 9mte). 7qvk was dropped
  at 3.1 Å, just above the cutoff.
- **Species bias**: most representatives come from llama (234) and alpaca (180),
  because these are the animals labs immunise. Models trained on this will see mostly
  camelid-derived sequences.
- **Clustering uses the full VH sequence.** Nanobody frameworks are very similar and
  binding is driven mostly by CDR-H3, so 13 of 551 clusters contain more than one
  CDR-H3 (mostly point mutants or library variants of the same binder).
- **The resolution filter removes ~44 % of rows**, many of them cryo-EM structures.
- **Snapshot**: SAbDab table downloaded September 2026 (latest entry June 2026).
- The raw SAbDab table (~17 MB) is not included; it has to be downloaded manually.
- The EDA notebook is committed without outputs.

## Reproducing

Requirements: Python ≥ 3.10 and `pip install -r requirements.txt`.

1. Download the full "all structures" summary table from
   [SAbDab](https://sabdab.opig.stats.ox.ac.uk) and save it as
   `data/raw/sabdab_nano_all.tsv` (comma- or tab-separated both work).
2. From the repository root:

```
python scripts/01_fetch_sabdab_nano.py      # select nanobodies bound to protein, flag HER2
python scripts/02_fetch_pdb_metadata.py     # method and resolution from RCSB (internet, ~5-10 min, resumable)
python scripts/03_quality_filter.py         # --max-res to change the resolution cutoff
python scripts/04_dedup_identity.py         # --identity to change the clustering threshold
python scripts/05_build_release.py --version v1.1 --note "what changed"
```

3. Open `notebooks/01_eda.ipynb` and run all cells.

Everything is pure Python (no CD-HIT or MMseqs2 install needed). Releases are never
overwritten; `05_build_release.py` refuses to reuse an existing version number.

## Repository layout

```
├── DECISIONS.md          curation log: every decision and why (Spanish)
├── CHANGELOG.md
├── scripts/              01–05, run in order
├── notebooks/01_eda.ipynb
└── data/
    ├── raw/pdb_metadata.csv    RCSB method/resolution per PDB entry (needed by the notebook)
    └── releases/v1.0/          the dataset
```

## Data sources

- **SAbDab / SAbDab-nano**, Oxford Protein Informatics Group.
  Dunbar J. et al. (2014) *SAbDab: the structural antibody database.* Nucleic Acids Res. 42:D1140–D1146.
  Schneider C., Raybould M.I.J., Deane C.M. (2022) *SAbDab in the age of biotherapeutics: updates including SAbDab-nano, the nanobody structure tracker.* Nucleic Acids Res. 50:D1368–D1372.
- **RCSB Protein Data Bank** Data API, for experimental method, resolution and deposition date.

If you use this dataset, please cite SAbDab and the PDB entries it is built from.
