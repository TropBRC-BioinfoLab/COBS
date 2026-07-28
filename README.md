# COBS

**COBS (Crisp to Overlap Boundary Similarity)** is a controlled overlapping-community detection framework. It starts from a greedy-modularity crisp partition, identifies boundary nodes, assigns additional memberships using node-to-community neighborhood similarity, and applies controlled community merging.

This repository contains the COBS implementation, a Streamlit interface, community-aware centrality analysis, comparison algorithms, and data used in the breast-cancer and CORUM experiments.

> **Project status:** research software / pre-release. The associated manuscript is in preparation and is intentionally not included in this public repository. Interfaces, file locations, and dependencies may still change before the first archived release.

## Main features

- Controlled crisp-to-overlap community refinement.
- Six neighborhood-similarity options: cosine/Salton, Jaccard, Sørensen–Dice, overlap coefficient, Adamic–Adar, and resource allocation.
- Overlapping community-consideration degree, closeness, and betweenness centralities.
- Structural evaluation, visualization, enrichment analysis, and hallmark-recovery analysis.
- Streamlit interface and a batch runner for PPI experiments.
- Reproducibility data for the breast-cancer and CORUM-derived networks.

## Repository structure

```text
.
├── app.py                         # Streamlit application
├── core/                          # Algorithms, metrics, centrality, data, and UI modules
│   └── algos/
│       └── cobs_overlap.py        # COBS implementation
├── scripts/
│   └── batch_run_overlap_ppi_v7_generic.py
├── data/
│   ├── breast_cancer/             # Breast-cancer network inputs and exports
│   └── corum/                     # CORUM-derived benchmark inputs
├── docs/
│   ├── DATA_SOURCES.md
│   └── REPRODUCIBILITY.md
├── CITATION.cff
├── CONTRIBUTING.md
├── CHANGELOG.md
└── requirements.txt
```

The Yeast benchmark files will be added after the original data, runner, parameters, and results have been consolidated and verified.

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
.venv\Scripts\Activate.ps1

# Linux/macOS
source .venv/bin/activate
```

Install the dependencies:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Some comparison algorithms distributed through CDlib may require additional optional packages. PyTorch is required by the included NOCD-GNN module because the current algorithm registry imports that module when the application starts.

## Run the Streamlit application

```bash
streamlit run app.py
```

The application opens with the Zachary Karate Club network. It can also build PPI networks from a disease query, a gene list, or a local identifier file by calling UniProt and STRING services. Internet access is required for those API-backed options.

For the included breast-cancer DisGeNET-derived identifier list, enter:

```text
data/breast_cancer/network/breastcancerdisgenet.txt
```

in the file-path field in the application.

## Run a PPI batch experiment

The current command-line workflow is the generic batch runner:

```bash
python scripts/batch_run_overlap_ppi_v7_generic.py \
  --input-file data/breast_cancer/network/breastcancerdisgenet.txt \
  --dataset-label breast_cancer \
  --species 9606 \
  --required-score 900 \
  --add-nodes 0 \
  --use-giant true \
  --include-cobs \
  --include-keys cobs \
  --timeout-sec 600 \
  --outdir outputs_breast_cancer
```

On Windows PowerShell, place the command on one line or replace each trailing `\` with a backtick.

Useful arguments include:

| Argument | Meaning | Default |
|---|---|---:|
| `--species` | NCBI taxonomy identifier | `9606` |
| `--required-score` | STRING confidence score, 0–999 | `400` |
| `--add-nodes` | Additional STRING interaction partners | `0` |
| `--use-giant` | Analyze the giant component (`auto`, `true`, `false`) | `auto` |
| `--include-cobs` | Include COBS in the batch | off |
| `--include-keys` | Comma-separated algorithm keys | all eligible |
| `--timeout-sec` | Per-algorithm timeout in seconds; `0` disables it | `300` |
| `--outdir` | Output directory | `outputs_breastcancer_batch` |

## COBS parameters

The COBS module exposes the following principal parameters:

| Parameter | Description | Default |
|---|---|---:|
| `sim_metric` | Neighborhood similarity: `cosine`, `jaccard`, `dice`, `overlap`, `aa`, or `ra` | `cosine` |
| `theta_sim` | Similarity threshold for overlap assignment and controlled merging | `0.40` |
| `theta_overlap` | Community overlap-ratio threshold used during merging | `0.30` |

The exact configuration for a reported experiment should be taken from its result files or execution manifest rather than inferred from these defaults.

## Included research data

- `data/breast_cancer/network/breast_cancer_network_export/breast_cancer_string_score900_gcc_edges.csv` is the final breast-cancer STRING GCC edge list used in the study (252 proteins and 729 unique interactions).
- `data/corum/corum_score900_gcc_5090_37942.csv` is the fixed CORUM-derived GCC input used by the benchmark (5,090 nodes and 37,942 edges).

These files are derived from third-party resources and do not automatically inherit the eventual software-code license. See [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md) before reusing or redistributing them.

## Reproducibility status

The repository currently provides the main code and the breast-cancer and CORUM inputs. It does not yet constitute the final archival reproduction package because the Yeast pipeline, final benchmark runner organization, environment lock file, and publication-linked outputs are still being consolidated. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). GitHub will expose them through **Cite this repository**. The journal article DOI will be added only after publication.

## License

A software license has not yet been selected by the project owners. Until a `LICENSE` file is added, copyright is retained and no open-source permission should be inferred. Third-party datasets remain subject to their respective source terms.

## Funding

This research was funded by the Indonesian Endowment Fund for Education (LPDP) on behalf of the Indonesian Ministry of Higher Education, Science and Technology, and managed under the EQUITY Program (Contract Nos. 4297/B3/DT.03.08/2025 and 42011/IT3/HK.07.00-4/P/B/2025).

## Contact

- TropBRC Bioinformatics Lab: <https://github.com/TropBRC-BioinfoLab>
- Corresponding researcher: Wisnu Ananta Kusuma, <ananta@apps.ipb.ac.id>

