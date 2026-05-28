# MEWS Usage Examples

This directory contains example scripts demonstrating how to use MEWS for equivalent width measurements.

## Available Examples

### `basic_usage.py` - Simple Single Galaxy Measurement

Demonstrates the core workflow:
- Loading spectrum data
- Selecting spectral lines
- Measuring equivalent widths
- Generating diagnostic plots

**Use this if:** You're new to MEWS and want to see how it works.

**Run it:**
```bash
python examples/basic_usage.py
```

---

## Batch Processing

For processing multiple galaxies, use the CLI interface:

### JSON Files (SDSS spectra)
```bash
python MEWS_Measuring_Equivalent_Width_in_Spectra.py \
  --preset example_json \
  --no-plots
```

### FITS Files (MaNGA data cubes)
```bash
python MEWS_Measuring_Equivalent_Width_in_Spectra.py \
  --input-type fits \
  --folder-path /path/to/fits/files \
  --redshift-csv /path/to/redshifts.csv \
  --output-csv results.csv
```

---

## Common Spectral Lines

### Balmer Series (Absorption in E+A galaxies)
- **H-alpha**: 6564.61 Å
- **H-beta**: 4862.68 Å
- **H-gamma**: 4341.68 Å
- **H-delta**: 4102.89 Å (strong E+A diagnostic)

### Emission Lines
- **[OIII]**: 5008.24 Å
- **[OII]**: 3729.875 Å
- **[NII]**: 6549.86, 6585.2 Å
- **[SII]**: 6718.295, 6732.674 Å

---

## Continuum Methods

MEWS supports 4 continuum fitting algorithms:

1. **`improved`** (default, recommended) - 5th-degree polynomial with outlier rejection
2. **`segmented`** - Fit continuum in segments with spline interpolation
3. **`robust`** - Iterative outlier rejection with MAD statistics
4. **`median`** - Median filtering approach

Select with: `--continuum-method improved`

---

## Output

MEWS generates:

1. **CSV file** with EW measurements for each galaxy
2. **Diagnostic plots** (optional) showing:
   - Continuum fit
   - Line profile
   - Integration region
   - Quality metrics

---

## Troubleshooting

**"Module not found" error:**
```bash
pip install -r requirements.txt
```

**No output plots:**
- Check that `--plot` flag is enabled (or `plot=True` in function call)
- Verify `save_dir` exists or can be created

**Poor continuum fits:**
- Try different continuum methods: `--continuum-method robust`
- Adjust width guess if line is very narrow/broad
- Check for contamination from nearby lines

---

## Questions?

See main README.md or contact: oliviaallegragreene@gmail.com

---

**Last Updated:** January 2026
