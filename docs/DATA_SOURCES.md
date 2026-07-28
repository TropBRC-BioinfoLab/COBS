# Data sources and reuse notes

The `data/` directory contains source-derived and processed research files. Their presence in this repository does not mean that they are covered by the eventual software-code license.

## Breast-cancer network

The breast-cancer analysis combines disease-associated identifiers with protein-interaction information. The repository currently includes:

- a disease-associated identifier list and related exported files;
- STRING identifier mappings and raw interaction responses;
- the processed STRING score-900 giant-connected-component edge list used in the study.

Primary resources:

- DisGeNET: <https://disgenet.com/>
- STRING: <https://string-db.org/>
- UniProt: <https://www.uniprot.org/>

Users must review and comply with the current terms, citation requirements, and redistribution rules of each source. In particular, raw or near-raw exports should not be assumed to be freely redistributable merely because they are present in a research workspace.

The processed study edge list is located at:

```text
data/breast_cancer/network/breast_cancer_network_export/
breast_cancer_string_score900_gcc_edges.csv
```

The validated study instance contains 252 unique proteins and 729 unique interactions, uses Homo sapiens (NCBI taxon 9606), STRING required score 900, no additional interactors, and only the giant connected component.

## CORUM-derived benchmark

CORUM is the Comprehensive Resource of Mammalian Protein Complexes:

- CORUM: <https://mips.helmholtz-muenchen.de/corum/>

The fixed benchmark input used in the reported run is:

```text
data/corum/corum_score900_gcc_5090_37942.csv
```

It contains 5,090 nodes and 37,942 edges. The original execution manifest should be retained unchanged as provenance even when its recorded input path points to the computer on which the run was performed. Repository documentation should map that historical path to the current relative location rather than altering the original record.

Check the CORUM download page for the applicable license and attribution requirements before redistribution or reuse. At the time of repository preparation, the project records CORUM as subject to CC BY-NC 4.0; users should verify that this remains current.

## Yeast benchmark

The Yeast benchmark data and its complete provenance are not yet included. They should be added only as a verified set comprising:

- the exact network input;
- the complex or group reference used for evaluation;
- the runner and parameters;
- the original manifest;
- the reported outputs;
- source and licensing information.

## Recommended citation practice

When publishing results produced with these data, cite both COBS and every underlying database or service used to construct the network. Database citations and access dates should be taken from the specific experiment metadata or manuscript, not inferred from this summary.

