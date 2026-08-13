# Changelog

All notable changes to this research software will be documented here. The format follows the principles of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers are intended to follow semantic versioning after the first stable release.

## [Unreleased]

### Changed

- Rebuilt the public data tree around fixed manuscript inputs and verified processed outputs.
- Replaced dataset-specific Streamlit file options with a generic user-supplied identifier-file workflow.
- Updated repository documentation for the final breast-cancer, CORUM-derived, and yeast study instances.
- Pinned the CYC2008 provenance link to the exact DyCluster commit used in the study.
- Made the generic batch runner import the repository package reliably when invoked through `scripts/`.
- Aligned identifier-file parsing with the documented line/comma/semicolon/tab-separated input format.
- Removed machine-local traceback artifacts from the public yeast benchmark package.
- Deferred software version and release-date metadata in `CITATION.cff` until a tagged publication release is created.

### Removed

- Removed the original DISGENET query export and derived disease-associated identifier lists from the public repository.
- Removed noncanonical CORUM raw/legacy files from the publication data tree.
- Removed the duplicated raw CYC2008 catalogue from the publication repository; the pinned upstream source is documented instead.

### Planned

- Consolidate the CORUM benchmark runner against the single canonical `core/` package.
- Add a focused COBS command-line entry point.
- Add automated tests and a reproducible environment lock file.
- Add a publication-linked release/tag and the article DOI after publication.
- Select and add the software license after owner approval.

## Initial public-repository preparation (July 2026)

### Added

- Initial public-repository preparation for COBS.
- Streamlit application and algorithm registry.
- COBS implementation and comparison algorithms.
- Overlapping community-consideration centrality modules.
- Generic PPI batch runner.
- Breast-cancer and CORUM-derived research inputs.
- Repository citation, contribution, data-source, and reproducibility documentation.

