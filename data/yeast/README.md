# Yeast benchmark dataset

This directory contains the yeast data used to evaluate COBS and the comparison
methods reported in the manuscript.

## Dataset provenance

- Organism: *Saccharomyces cerevisiae*
- STRING taxonomy ID: `4932`
- Reference catalogue: CYC2008
- Catalogue size: 408 manually curated protein complexes
- Scientific reference: Pu S, Wong J, Turner B, Cho E, Wodak SJ.
  "Up-to-date catalogues of yeast protein complexes." *Nucleic Acids Research*
  2009;37(3):825-831. https://doi.org/10.1093/nar/gkn1005
- Archived catalogue copy:
  https://github.com/emhanna/DyCluster/blob/master/Supplementary%20Materials/Datasets/CYC2008%20Catalogue.txt
- Archive accessed: 2026-07-28

The original CYC2008 website was no longer accessible. The catalogue file in
this directory is therefore a publicly archived copy from the DyCluster
repository. It contains the same 408 complexes described in the original paper.

## Directory structure

- `CYC2008 Catalogue.txt`: original catalogue copy containing complexes C1-C408.
- `cyc2008_unique_proteins.txt`: 1,627 unique protein identifiers extracted
  from the catalogue.
- `string_400/`, `string_700/`, and `string_900/`: STRING-derived networks,
  identifier mappings, node attributes, run metadata, full COBS parameter
  searches, top-20 parameter results, and best COBS memberships for each
  confidence threshold.
- `benchmark_score900/`: comparative benchmark summary for the final
  0.9-confidence analysis, together with the COBS communities, memberships,
  run metadata, and error log.
- `SHA256SUMS`: SHA-256 checksums for integrity verification.

## Final configuration reported in the manuscript

- STRING required score: 900
- Giant component: 1,531 nodes and 26,721 edges
- COBS similarity metric: Adamic-Adar (`AA`)
- Similarity threshold: 0.5
- Overlap threshold: 0.1
- Detected communities: 29
- Shen modularity: 0.5821
- Conductance: 0.1547

The score-900 `edges.csv` contains 26,844 retrieved interactions. Of these,
26,721 connect pairs within the 1,531-node giant component used as the active
benchmark graph; 123 interactions involve nodes outside that component.

The catalogue complexes define the initial curated yeast protein set. The
community-detection methods were evaluated on the STRING-derived interaction
network; the 408 catalogue complexes are not COBS output communities.
