# Reproducibility guide

This document distinguishes what can currently be reproduced from this repository from components that are still being consolidated.

## Current coverage

Available:

- COBS implementation and default parameters;
- Streamlit analysis interface;
- generic PPI batch runner;
- breast-cancer disease identifier and network-export files;
- final breast-cancer score-900 GCC edge list (252 nodes, 729 edges);
- fixed CORUM-derived GCC input (5,090 nodes, 37,942 edges);
- community, metric, centrality, visualization, enrichment, and hallmark-analysis modules.

Pending:

- complete Yeast network, reference data, runner, parameters, manifest, and outputs;
- a consolidated CORUM benchmark entry point in this public folder;
- publication-linked final result tables;
- a fully locked environment and automated regression tests.

## Preserve original provenance

Execution manifests are historical records. Do not rewrite machine-specific paths, timestamps, platform details, or software versions in an original manifest. If a file has moved, document the mapping separately. For the current CORUM package, the canonical repository copy of the fixed input is:

```text
data/corum/corum_score900_gcc_5090_37942.csv
```

## Minimum information for a reproducible run

Record the following for every experiment:

1. input filename and SHA-256 checksum;
2. node and edge counts before and after component filtering;
3. species/taxon and data-source versions or retrieval dates;
4. STRING required score and `add_nodes` setting, if used;
5. graph scope (full network or giant component);
6. algorithm key and all parameter values;
7. random seed, if supported;
8. timeout and post-processing settings;
9. Python, dependency, and operating-system versions;
10. output files, runtime, status, and timestamp.

## Breast-cancer study instance

The final edge list can be validated independently with the following expectations:

| Property | Expected value |
|---|---:|
| Unique nodes | 252 |
| Unique undirected edges | 729 |
| Self-loops | 0 |
| Duplicate undirected edges | 0 |
| Species | Homo sapiens (taxon 9606) |
| STRING required score | 900 |
| Additional interactors | 0 |
| Graph scope | Giant connected component |

The manuscript's selected COBS configuration and all reported result values must be read from the final experiment outputs, not inferred from software defaults.

## Environment capture

After installing a clean environment and successfully running the relevant workflows, capture the resolved environment without replacing `requirements.txt`:

```bash
python --version
python -m pip freeze > environment-lock.txt
```

The lock file should be created from the environment actually used for the archived release. Do not generate it from an unrelated development environment.

