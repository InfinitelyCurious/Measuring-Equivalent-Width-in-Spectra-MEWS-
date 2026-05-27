# MEWS: Measuring Equivalent Width in Spectra

An automated pipeline for accurate equivalent width measurements in galaxy spectra, with built-in validation and adaptive continuum fitting algorithms.

## Overview

MEWS addresses systematic measurement errors in automated spectral analysis pipelines by combining multiple continuum-fitting algorithms with visual verification capabilities. Originally developed for post-starburst galaxy classification in SDSS-IV MaNGA, MEWS is generalizable to any spectroscopic dataset requiring robust equivalent width measurements.

## Key Features

- **Adaptive Continuum Fitting**: Multiple algorithms (5th-degree polynomial, segmented, robust, median filter) automatically selected based on spectral characteristics
- **Smart Profile Fitting**: Automatic selection between Gaussian and Voigt profiles based on line depth (>20% continuum uses Voigt)
- **Built-in Validation**: Visual diagnostic plots at each measurement step
- **Flexible Input**: Handles both `.json` (SDSS single-fiber) and `.fits` (MaNGA IFS) formats
- **Automated Quality Control**: Detects unrealistic fits, out-of-range lines, and failed measurements
- **Wavelength-Adaptive Windows**: 
  - 100 Å for red lines (≥6500 Å): Hα, [S II], [N II]
  - 80 Å for mid-spectrum (≥4800 Å): Hβ, [O III]
  - 50 Å for blue lines (<4800 Å): Hγ, Hδ, [O II]

## Why MEWS?

Standard spectral synthesis models systematically **underestimate** deep Balmer absorption lines characteristic of post-starburst galaxies. MEWS was developed to:

1. Reproduce the accuracy of manual measurements at scale
2. Provide transparent, auditable continuum fitting decisions
3. Enable systematic equivalent width analysis across large samples (100s-1000s of galaxies)

**Validation**: MEWS successfully reproduced manual PyRAF measurements for 30 E+A galaxies (Greene et al. 2021), then scaled to measure 579+ candidates in MPL-11.

## Installation

```bash
# Clone repository
git clone https://github.com/InfinitelyCurious/Measuring-Equivalent-Width-in-Spectra-MEWS-.git
cd Measuring-Equivalent-Width-in-Spectra-MEWS-

# Install dependencies
pip install -r requirements.txt
