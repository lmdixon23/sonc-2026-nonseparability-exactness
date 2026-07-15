# Tested release environment

## Python

- Python 3.13.5
- SymPy 1.14.0
- python-flint 0.9.0
- NumPy 2.4.4
- SciPy 1.17.1
- mpmath 1.3.0

Install with:

```bash
python -m pip install -r requirements-lock.txt
```

## TeX and PDF tools

The tested build used:

- pdfTeX 3.141592653-2.6-1.40.26
- TeX Live 2025/dev
- `pdfinfo`
- `pdffonts`
- `pdftotext`

The manuscript uses the packages `geometry`, `amsmath`, `amssymb`, `amsthm`, `mathtools`, `microtype`, `enumitem`, `booktabs`, `longtable`, `array`, `hyperref`, and `xcolor`.
