# Contributing

COBS is currently research software under active consolidation. Contributions that improve correctness, reproducibility, documentation, testing, and portability are welcome after coordination with the maintainers.

## Before opening a change

1. Open an issue describing the problem, proposed change, affected dataset or algorithm, and expected behavior.
2. Keep scientific-method changes separate from refactoring or interface changes.
3. Do not add manuscript drafts, credentials, private paths, restricted data, generated caches, or large intermediate outputs.
4. Confirm that any new dataset may legally be redistributed and document its source, version, retrieval date, transformation, and license or terms.

## Development setup

```bash
python -m venv .venv
pip install -r requirements.txt
```

Before submitting a change, verify at minimum:

```bash
python -m compileall -q app.py core scripts
```

Also run the affected workflow on a small graph and report:

- Python and operating-system versions;
- exact command or application settings;
- random seed, if applicable;
- input checksum or unambiguous data version;
- observed output and expected output.

## Coding and scientific-integrity expectations

- Prefer deterministic behavior where the underlying method permits it.
- Preserve original execution manifests; add explanatory metadata instead of rewriting provenance.
- Use relative repository paths in new scripts and documentation.
- Document parameter defaults and any changes that could alter published metrics.
- Add or update validation checks when changing graph construction, membership handling, metrics, or centrality calculations.
- Avoid committing `__pycache__`, `.pyc`, virtual environments, API credentials, or unreviewed third-party exports.

## Pull requests

A pull request should contain a concise rationale, files changed, validation performed, and any effect on previously reported results. Scientific changes should include before/after metrics and an explanation of why the new behavior is correct.

