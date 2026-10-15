# Computational companion workflow

- Keep verification runnable offline using the standard library. The optional source PDF check may require pypdf.
- Run `python3 scripts/check_reproduction.py` after changing the payload.
- Do not use Python optimization; assertions are part of verification.
- Preserve source citations and distinguish exact matching from numerical checks.
- Update manifest hashes deliberately when payloads change; retain source provenance.
- No manuscript PDF generation hook applies to this independent repository.
- Do not invent a DOI, remote URL, release version, license grant, or confirmed software authorship.
- Keep manuscript sources, figure assets and figure-generation scripts in the private manuscript repository. Never add them here.
- Manuscript comparison must remain optional via --manuscript.
