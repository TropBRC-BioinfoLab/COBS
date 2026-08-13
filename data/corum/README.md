# CORUM-derived structural benchmark

This directory contains the fixed processed network used for the final CORUM-derived structural benchmark.

## Canonical benchmark input

`corum_score900_gcc_5090_37942.csv`

- graph scope: giant connected component;
- STRING required score: 900;
- projected STRING species background: *Homo sapiens* (`NCBI taxon 9606`);
- nodes: 5,090;
- edges: 37,942.

This is the canonical processed input used by the final fixed benchmark. Earlier exploratory full-graph and noncanonical CORUM-derived files are intentionally not included in the publication repository.

## Source provenance

The benchmark was derived from the **Complete complexes** dataset of the **Comprehensive Resource of Mammalian Protein Complexes (CORUM)**, release **5.3**. The raw Complete complexes file is not redistributed in this repository.

Official CORUM download repository:

<https://mips.helmholtz-muenchen.de/corum/download/>

Direct current Complete-complexes TXT download endpoint recorded for this study:

<https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=complete&file_format=txt>

The Complete complexes collection contains annotated mammalian complexes from multiple organisms. CORUM-derived identifiers were mapped onto a common human STRING PPI background for structural benchmarking. The benchmark should therefore be interpreted as a **mammalian-derived structural reference projected onto the human STRING interactome**, not as a direct cross-species validation.

CORUM states that copyrightable database content is made available under **CC BY-NC 4.0**, subject to the additional disclaimer and rights conditions provided on the official site.
