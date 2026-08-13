# Data sources and reuse notes

The `data/` directory contains processed research files retained to support traceability and reproduction of the manuscript analyses. Their presence in this repository does **not** mean that they are covered by any eventual software-code license.

## Breast-cancer network

### Disease-associated source

Breast-cancer-associated genes were obtained from **DISGENET v26.1**, queried on **1 May 2026** for *Malignant neoplasm of breast* (`C0006142`) using the GDA view, curated source, no additional filtering, and no custom score. The query returned 755 records.

- DISGENET: <https://disgenet.com/>
- Repository provenance record: `data/breast_cancer/search_info_C0006142.txt`

The original query export and the extracted disease-associated identifier list are intentionally **not redistributed**. DISGENET states that its standard licenses do not permit redistribution or resale of the entire database or a portion of it. Users who need the original disease-associated data should obtain it directly from DISGENET under the applicable license.

### Protein associations

The final PPI network was obtained from **STRING v12.0** for *Homo sapiens* (`NCBI taxon 9606`) with `required_score = 900`, `add_nodes = 0`, followed by retention of the giant connected component.

- STRING v12.0: <https://version-12-0.string-db.org/>
- STRING information/licensing: <https://www.string-db.org/cgi/info>

STRING states that its data are available under **CC BY 4.0**, with attribution required.

The authoritative fixed manuscript input is:

```text
data/breast_cancer/breast_cancer_string_score900_gcc_edges.csv
```

It contains 252 proteins and 729 unique undirected interactions.

## CORUM-derived benchmark

The mammalian structural reference was derived from the **Complete complexes** dataset of the **Comprehensive Resource of Mammalian Protein Complexes (CORUM)**, release **5.3**.

- CORUM official site: <https://mips.helmholtz-muenchen.de/corum/>
- CORUM download repository: <https://mips.helmholtz-muenchen.de/corum/download/>
- Direct current Complete-complexes TXT endpoint recorded for this study: <https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=complete&file_format=txt>

The raw Complete complexes file is not redistributed. The final processed benchmark input is:

```text
data/corum/corum_score900_gcc_5090_37942.csv
```

It contains 5,090 nodes and 37,942 edges and is the fixed score-900 GCC used in the reported benchmark. CORUM-derived identifiers were projected onto a common human STRING background (`NCBI taxon 9606`), so this is a mammalian-derived structural reference on a human STRING interactome rather than a direct cross-species evaluation.

CORUM states that copyrightable database content is made available under **CC BY-NC 4.0**, subject to the additional disclaimer and rights conditions on the official download site.

## Yeast benchmark

The yeast reference set was based on the **CYC2008** catalogue of 408 manually curated *Saccharomyces cerevisiae* protein complexes described by Pu et al.:

- Publication: <https://doi.org/10.1093/nar/gkn1005>
- Exact archived catalogue used in the study, pinned to DyCluster commit `6d0fd9e1dc26a6b13dbf266c05a2c88fcc0c7f65`: <https://github.com/emhanna/DyCluster/blob/6d0fd9e1dc26a6b13dbf266c05a2c88fcc0c7f65/Supplementary%20Materials/Datasets/CYC2008%20Catalogue.txt>

The raw catalogue is not duplicated here. Study preprocessing produced 1,627 unique yeast identifiers, which were mapped to **STRING v12.0** using `NCBI taxon 4932`. Processed score-400, score-700, and score-900 networks and the final score-900 benchmark package are retained under `data/yeast/`.

## Other reference resources used in the manuscript

These resources are not redistributed in this repository:

- Cancer Hallmarks Core Cancer Hallmarks Gene Set: <https://cancerhallmarks.com/download_file/Menyhart_JPA_CancerHallmarks_core.txt>
- Cancer Hallmarks reference publication: <https://doi.org/10.1016/j.jpha.2024.101065>
- Enrichr: <https://maayanlab.cloud/Enrichr/>
- Gene Ontology: <https://geneontology.org/>
- KEGG: <https://www.kegg.jp/>
- UniProt: <https://www.uniprot.org/>

## Reuse principle

When reusing repository data, cite both COBS and the underlying databases/services from which each processed file was derived. Source terms, licenses, versions, and access rules can change; users are responsible for checking the current conditions at the original source.
