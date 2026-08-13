# Yeast benchmark data

This directory contains the processed yeast data used to evaluate COBS and comparison methods reported in the manuscript.

## Dataset provenance

- Organism: *Saccharomyces cerevisiae*
- STRING taxonomy ID: `4932`
- Reference catalogue: CYC2008
- Catalogue size: 408 manually curated protein complexes
- Scientific reference: Pu S, Wong J, Turner B, Cho E, Wodak SJ. “Up-to-date catalogues of yeast protein complexes.” *Nucleic Acids Research* 2009;37(3):825–831. <https://doi.org/10.1093/nar/gkn1005>

Because the original CYC2008 website is no longer accessible, the study used the public archived catalogue in the DyCluster repository, pinned to Git commit `6d0fd9e1dc26a6b13dbf266c05a2c88fcc0c7f65`:

<https://github.com/emhanna/DyCluster/blob/6d0fd9e1dc26a6b13dbf266c05a2c88fcc0c7f65/Supplementary%20Materials/Datasets/CYC2008%20Catalogue.txt>

The raw archived catalogue is not duplicated in this publication repository. The study preprocessing extracted **1,627 unique protein identifiers**, preserved here as `cyc2008_unique_proteins.txt`, and mapped them to STRING taxonomy ID 4932.

## Directory structure

- `cyc2008_unique_proteins.txt` — 1,627 unique protein identifiers extracted during study preprocessing.
- `string_400/`, `string_700/`, `string_900/` — processed STRING networks, identifier mappings, node attributes, run metadata, full COBS parameter searches, top-20 parameter results, and best COBS memberships for each confidence threshold.
- `benchmark_score900/` — comparative benchmark summary for the final score-900 analysis, including COBS communities, memberships, run metadata, summary metrics, and workbook.
- `SHA256SUMS` — SHA-256 checksums for the processed files in this directory.

The original comparative-run metadata is preserved as historical provenance. It may refer to an `errors.csv` artifact produced during execution. That traceback artifact is intentionally not redistributed in the publication repository because it contains machine-local filesystem paths; algorithm status and reported benchmark metrics remain available in `benchmark_score900/summary_metrics.csv` and `yeast_benchmark_results.xlsx`.

## Networks and selected configurations

| STRING required score | Active GCC | Selected similarity | `theta_sim` | `theta_overlap` | Shen modularity |
|---:|---:|---|---:|---:|---:|
| 400 | 1,627 proteins / 69,357 interactions | Adamic–Adar | 0.5 | 0.1 | 0.4188 |
| 700 | 1,613 proteins / 41,533 interactions | Adamic–Adar | 0.4 | 0.1 | 0.4998 |
| 900 | 1,531 proteins / 26,721 interactions | Adamic–Adar | 0.5 | 0.1 | 0.5821 |

For score 900, `string_900/edges.csv` contains 26,844 retrieved interactions before active-component filtering. The manuscript benchmark uses the 1,531-node giant component containing 26,721 interactions.

## Grid-search versus comparative-benchmark metrics

The score-specific `run_metadata.json` files record metrics generated during the COBS parameter-search workflow. `benchmark_score900/summary_metrics.csv` records the separate comparative benchmark execution used for the manuscript comparison table. The same selected score-900 COBS configuration and Shen modularity are reproduced in both workflows, while some secondary metrics can differ because they were produced by different evaluation stages. For manuscript comparative values, use `benchmark_score900/summary_metrics.csv` and `yeast_benchmark_results.xlsx`.
