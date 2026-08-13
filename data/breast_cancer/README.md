# Breast-cancer study data

This directory contains the public, processed data retained for the final breast-cancer PPI analysis.

## Files

- `breast_cancer_string_score900_gcc_edges.csv` — exact undirected STRING-derived giant-connected-component edge list used for the final analysis; 252 proteins and 729 unique interactions.
- `breast_cancer_string_score900_gcc_edges_metadata.json` — validation metadata and SHA-256 checksum for the final edge list.
- `search_info_C0006142.txt` — query-provenance record for the disease selection performed in DISGENET.
- `SHA256SUMS` — checksums for the public data files in this directory.

## Disease-query provenance

DISGENET was queried on 1 May 2026 for **Malignant neoplasm of breast** (UMLS CUI `C0006142`) using the Diseases analysis, GDA tab, curated source, no additional filter, and no custom score. The query returned **755 records**.

The original DISGENET query workbook and the extracted disease-associated identifier list are **not redistributed** in this repository. DISGENET data remain subject to DISGENET licensing conditions. The query-selection record is retained so the study provenance remains auditable without republishing the licensed query output.

DISGENET: <https://disgenet.com/>

## Final STRING network

The final reported breast-cancer network used:

- species: *Homo sapiens* (`NCBI taxon 9606`);
- STRING version: 12.0;
- `required_score = 900`;
- `add_nodes = 0`;
- giant connected component only;
- 252 proteins and 729 unique undirected interactions.

Version-specific STRING resource: <https://version-12-0.string-db.org/>

The fixed edge list in this directory is the authoritative network input for reproducing the final manuscript analysis independently of future changes to external APIs.
