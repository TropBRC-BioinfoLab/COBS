# COBS

**COBS (Crisp-to-Overlap Boundary Similarity)** is a controlled overlapping-community detection framework for protein-protein interaction (PPI) networks. It starts from a greedy-modularity crisp partition, identifies boundary nodes, assigns additional memberships using node-to-community neighborhood similarity, and applies controlled community merging.

This repository is the public research-software and processed-data companion to the manuscript **“COBS: controlled overlapping community detection and overlapping community-consideration centrality for breast cancer protein–protein interaction networks.”**

## What is included

- the COBS implementation and comparison-algorithm registry;
- overlapping community-consideration centrality modules;
- a Streamlit analysis interface;
- a generic command-line PPI batch runner;
- the exact processed breast-cancer STRING score-900 GCC used for the final case study;
- the fixed CORUM-derived score-900 GCC used for the final structural benchmark;
- processed yeast STRING networks, COBS parameter-search outputs, and the final score-900 benchmark package;
- data-source and reproducibility documentation.

Raw third-party database exports that are not required for redistribution are intentionally **not** bundled. In particular, the original DISGENET query export, the CORUM Complete complexes source file, and the archived CYC2008 catalogue are referenced at their source locations rather than redistributed here.

## Repository structure

```text
.
├── app.py
├── core/
│   └── algos/
│       └── cobs_overlap.py
├── scripts/
│   └── batch_run_overlap_ppi_v7_generic.py
├── data/
│   ├── breast_cancer/
│   ├── corum/
│   └── yeast/
├── docs/
│   ├── DATA_SOURCES.md
│   ├── PUBLICATION_DATA_MAP.md
│   └── REPRODUCIBILITY.md
├── CITATION.cff
├── CONTRIBUTING.md
├── CHANGELOG.md
└── requirements.txt
```

## Installation

Python 3.11 or 3.12 is recommended.

```bash
git clone https://github.com/TropBRC-BioinfoLab/COBS.git
cd COBS
python -m venv .venv
```

Activate the environment:

```bash
# Windows PowerShell
.\.venv\Scripts\Activate.ps1

# Linux/macOS
source .venv/bin/activate
```

Install the dependencies:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Some comparison algorithms exposed through CDlib can require additional optional packages. PyTorch is required by the included NOCD-GNN module because the current algorithm registry imports that module when the application starts.

## Run the Streamlit application

```bash
streamlit run app.py
```

The application opens with the Zachary Karate Club graph and can also construct PPI networks from:

- a disease query through UniProt and STRING;
- a pasted gene/protein identifier list;
- a **user-supplied local identifier file**.

The repository does not bundle licensed disease-query exports. If an analysis starts from a licensed third-party identifier set, obtain it from the original provider under the applicable terms and keep it outside the public repository unless redistribution is permitted.

## Run a generic PPI batch experiment

The command-line runner builds a STRING network from a user-supplied identifier file:

```bash
python scripts/batch_run_overlap_ppi_v7_generic.py \
  --input-file identifiers.txt \
  --dataset-label my_dataset \
  --species 9606 \
  --required-score 900 \
  --add-nodes 0 \
  --use-giant true \
  --include-cobs \
  --include-keys cobs \
  --timeout-sec 600 \
  --outdir outputs_ppi_batch
```

On Windows PowerShell, place the command on one line or replace each trailing `\` with a backtick.

Useful arguments include:

| Argument | Meaning | Default |
|---|---|---:|
| `--input-file` | Local identifier file | `identifiers.txt` |
| `--species` | NCBI taxonomy identifier | `9606` |
| `--required-score` | STRING confidence score, 0–999 | `400` |
| `--add-nodes` | Additional STRING interaction partners | `0` |
| `--use-giant` | Analyze the giant component (`auto`, `true`, `false`) | `auto` |
| `--include-cobs` | Include COBS in the batch | off |
| `--include-keys` | Comma-separated algorithm keys | all eligible |
| `--timeout-sec` | Per-algorithm timeout in seconds; `0` disables it | `300` |
| `--outdir` | Output directory | `outputs_ppi_batch` |

External database services can change over time. For reproduction of the reported manuscript analyses, use the **fixed processed networks and archived outputs** supplied here and in Supplementary Data S1 rather than assuming a future API query will recreate byte-identical source data.

## COBS parameters

The COBS module exposes the following principal parameters:

| Parameter | Description | Default |
|---|---|---:|
| `sim_metric` | Neighborhood similarity: `cosine`, `jaccard`, `dice`, `overlap`, `aa`, or `ra` | `cosine` |
| `theta_sim` | Similarity threshold for overlap assignment and controlled merging | `0.40` |
| `theta_overlap` | Community overlap-ratio threshold used during merging | `0.30` |

Software defaults are not the reported manuscript settings. The exact configurations used for reported experiments are recorded in the processed data and publication map.

## Manuscript study instances

| Study setting | Final graph | Selected COBS configuration | Reported Shen modularity |
|---|---:|---|---:|
| CORUM-derived score 900 | 5,090 nodes / 37,942 edges | cosine, `theta_sim=0.3`, `theta_overlap=0.1` | 0.7026 |
| Yeast score 900 | 1,531 nodes / 26,721 edges | Adamic–Adar, `theta_sim=0.5`, `theta_overlap=0.1` | 0.5821 |
| Breast cancer score 900 | 252 nodes / 729 edges | Adamic–Adar, `theta_sim=0.3`, `theta_overlap=0.1` | 0.5909 |

See [`docs/PUBLICATION_DATA_MAP.md`](docs/PUBLICATION_DATA_MAP.md) for the relationship between repository files, Supplementary Data S1, and manuscript results.

## Data licensing and provenance

Processed study files in `data/` are derived from third-party resources and do not automatically inherit the software-code licensing status of this repository. See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) before reusing or redistributing them.

## Reproducibility status

The repository provides the main research code and the principal processed inputs used in the reported analyses. The exact breast-cancer final edge list, fixed CORUM benchmark input, and processed yeast benchmark package are included. A fully locked environment, automated regression tests, and a dedicated CORUM benchmark entry point remain future repository-hardening tasks. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). The article DOI can be added after publication.

## License

A software license has not yet been selected by the project owners. Until a `LICENSE` file is added, copyright is retained and no open-source permission should be inferred. Third-party data remain subject to their respective source terms and licenses.

## Funding

This research was funded by the Indonesian Endowment Fund for Education (LPDP) on behalf of the Indonesian Ministry of Higher Education, Science and Technology, and managed under the EQUITY Program (Contract Nos. 4297/B3/DT.03.08/2025 and 42011/IT3/HK.07.00-4/P/B/2025).

## Contact

- TropBRC Bioinformatics Lab: <https://github.com/TropBRC-BioinfoLab>
- Corresponding researcher: Wisnu Ananta Kusuma, <ananta@apps.ipb.ac.id>
