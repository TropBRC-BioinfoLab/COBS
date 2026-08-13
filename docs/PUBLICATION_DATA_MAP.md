# Publication data map

This file maps the principal manuscript study settings to the public repository and Supplementary Data S1. It is intended to prevent confusion between exploratory files and the fixed inputs/results used in the final manuscript.

## Summary

| Study setting | Fixed graph used in manuscript | Selected COBS parameters | Repository | Supplementary Data S1 |
|---|---:|---|---|---|
| CORUM-derived score 900 | 5,090 nodes / 37,942 edges | cosine; `theta_sim=0.3`; `theta_overlap=0.1` | `data/corum/corum_score900_gcc_5090_37942.csv` | `corum_summary_metrics_fixed.csv`, `corum_batch_results_fixed.xlsx` |
| Yeast score 900 | 1,531 nodes / 26,721 edges | AA; `theta_sim=0.5`; `theta_overlap=0.1` | `data/yeast/string_900/`, `data/yeast/benchmark_score900/` | `yeast_benchmark_results.xlsx`, score-900 yeast grid-search workbook |
| Breast cancer score 900 | 252 nodes / 729 edges | AA; `theta_sim=0.3`; `theta_overlap=0.1` | `data/breast_cancer/breast_cancer_string_score900_gcc_edges.csv` | exact edge list, metadata, grid-search output, memberships, centrality, enrichment, hallmark outputs |

## CORUM-derived benchmark

Raw source: **CORUM Complete complexes**, release 5.3. The raw file is external-only. The publication repository contains only the canonical fixed score-900 GCC used in the final benchmark.

Reported COBS result: Shen modularity **0.7026**, 59 communities, conductance **0.1915**, average transitivity **0.3962** in the comparative benchmark.

Do not substitute earlier exploratory full-graph CORUM outputs or `humanComplexes` data for the fixed manuscript benchmark.

## Yeast benchmark

Raw reference: CYC2008, 408 manually curated complexes; study preprocessing extracted 1,627 unique identifiers. The raw archived catalogue remains at the pinned DyCluster source and is not duplicated in this repository.

Processed STRING networks:

| Required score | Active GCC | Selected COBS parameters | Shen modularity |
|---:|---:|---|---:|
| 400 | 1,627 / 69,357 | AA, 0.5 / 0.1 | 0.4188 |
| 700 | 1,613 / 41,533 | AA, 0.4 / 0.1 | 0.4998 |
| 900 | 1,531 / 26,721 | AA, 0.5 / 0.1 | 0.5821 |

For the manuscript comparative table at score 900, use `data/yeast/benchmark_score900/summary_metrics.csv` and `yeast_benchmark_results.xlsx`; the reported COBS conductance is **0.1547**.

## Breast-cancer case study

Disease-source query: DISGENET v26.1, queried 1 May 2026 for *Malignant neoplasm of breast* (`C0006142`), GDA view, curated source, no additional filter/custom score; **755 records** returned.

The original DISGENET query export and extracted identifier list are not redistributed. The public repository retains only `search_info_C0006142.txt` as provenance and the fixed processed STRING network.

STRING sensitivity graphs reported in the manuscript:

| Required score | GCC |
|---:|---:|
| 400 | 651 proteins / 6,279 interactions |
| 700 | 433 proteins / 2,034 interactions |
| 900 | 252 proteins / 729 interactions |

The final case study uses score 900. Selected COBS parameters: Adamic–Adar, `theta_sim = 0.3`, `theta_overlap = 0.1`. Reported result: Shen modularity **0.5909**, 12 communities, with BRCA1 and STK11 assigned to two communities.

The complete final memberships, OCC centrality, enrichment outputs, hallmark-recovery results, and breast-cancer grid-search workbooks are archived in Supplementary Data S1 rather than duplicated throughout the code repository.

## Source-of-truth rule

If similarly named historical files disagree, use this order:

1. fixed publication input/output identified in this map and Supplementary Data S1;
2. final manuscript values;
3. repository run metadata for the corresponding fixed workflow;
4. exploratory/legacy files only for historical investigation, not for publication reporting.
