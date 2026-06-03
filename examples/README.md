# MEWS Usage Examples

This directory contains example scripts demonstrating how to use MEWS for equivalent width measurements.

## Available Examples

### 1. Basic Usage (`basic_usage.py`)

Simple example that:
- Loads data from configured paths
- Applies dust correction
- Creates a color-mass diagram

**Use this if:** You want the quickest way to generate a plot with default settings.

### 2. Custom Sample Analysis (`custom_sample_analysis.py`)

Advanced example that:
- Loads and inspects data
- Analyzes green valley membership
- Demonstrates customization options

**Use this if:** You want to understand the workflow and customize analysis.

---

## Running Examples

### Prerequisites

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Download data files** (see `data/README.md`)

3. **Update configuration** in your Python script or use command-line interface

### Run Basic Example

```bash
python examples/basic_usage.py
```

### Run Custom Analysis

```bash
python examples/custom_sample_analysis.py
```

---

## Expected Output

Both examples will:
1. Print status messages to console
2. Generate plot(s) in `output/` directory (if plotting enabled)
3. Display summary statistics

---

## Diagnostic Plots & Quality Control

MEWS generates detailed visual diagnostics for every measurement, enabling quality assessment at scale.

### What You'll See:

For each galaxy, MEWS can produce:

1. **Full spectrum continuum fit** - Shows overall continuum placement across entire wavelength range
2. **Per-line diagnostics** (for each of 10 spectral lines):
   - Local continuum fit around the target line
   - Line profile fit (Gaussian or Voigt) with integration region
   - Normalized spectrum showing line depth/strength
   - Residuals to identify poor fits

**Total diagnostic output**: Up to 40+ plots per galaxy (1 full spectrum + ~4 plots × 10 lines)

### How MEWS Was Used for Catalog Curation:

```python
# Process candidate sample with full diagnostics
for galaxy in candidate_list:
    measure_ew(
        file_path=f"data/{galaxy['plateifu']}.fits",
        file_type='fits',
        line_name='Hdelta',
        plot=True,  # Enable diagnostic plots
        save_dir=f"diagnostics/{galaxy['plateifu']}/",
        continuum_method='improved'
    )
```

**Practical application**: Visual inspection of these diagnostic plots enabled curation of ~600 E+A galaxy candidates → 183 validated systems with high-quality measurements. Each galaxy's diagnostic suite revealed:
- Continuum fitting quality
- Line profile appropriateness (Gaussian vs. Voigt)
- Integration boundary placement
- Spectral artifacts or contamination

This visual validation was essential for maintaining measurement quality comparable to manual inspection while processing hundreds of galaxies.

### Disabling Plots for Production Runs:

```python
# Fast processing without diagnostic output
measure_ew(
    ...,
    plot=False  # Disable plots for speed
)
```

**Note**: MEWS was specifically designed to *replace* obsolete PyRAF tools. Initial validation against 30 manual PyRAF measurements (Greene et al. 2021) showed excellent agreement, but PyRAF is no longer maintained. MEWS provides modern, reproducible equivalent width measurements with superior diagnostic capabilities.

---

## Troubleshooting

**"Module not found" error:**
```bash
pip install -r requirements.txt
```

**No output plots:**
- Check that `plot=True` in function call
- Verify `save_dir` path exists or can be created
- Check disk space (diagnostic plots can accumulate quickly!)

**Poor continuum fits:**
- Try different continuum methods: `continuum_method='robust'`
- Adjust `width_guess` if line is very narrow/broad
- Review diagnostic plots to identify contamination from nearby lines
- Check for artifacts in spectrum (cosmic rays, bad pixels)

**"Line out of wavelength range" warnings:**
- Some SDSS/MaNGA spectra have truncated blue or red ends
- Check wavelength coverage before selecting lines to measure
- MEWS will skip lines outside spectral range and note in output

---

## Questions?

See main README.md or contact: oliviaallegragreene@gmail.com

---

**Last Updated:** January 2026
```
