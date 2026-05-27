MEWS: Measuring Equivalent Width in Spectra

An automated pipeline for accurate equivalent width measurements in galaxy spectra, with built-in validation and adaptive continuum fitting algorithms.

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

MEWS addresses systematic measurement errors in automated spectral analysis pipelines by combining multiple continuum-fitting algorithms with visual verification capabilities. Originally developed for post-starburst galaxy classification in SDSS-IV MaNGA, MEWS is generalizable to any spectroscopic dataset requiring robust equivalent width measurements.

Developed as part of dissertation research at Vanderbilt University (2020-2026), MEWS successfully processed 579+ galaxy spectra for the largest spatially-resolved E+A galaxy catalog to date.

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

**Validation**: MEWS successfully reproduced manual PyRAF measurements for 30 E+A galaxies (Greene et al. 2021, ApJ), then scaled to measure 579+ candidates in MaNGA MPL-11.

## Installation

```bash
# Clone repository
git clone https://github.com/InfinitelyCurious/Measuring-Equivalent-Width-in-Spectra-MEWS-.git
cd Measuring-Equivalent-Width-in-Spectra-MEWS-

# Install dependencies
pip install -r requirements.txt
Dependencies
astropy>=5.0 (spectral models, sigma clipping)
specutils>=1.0 (Spectrum1D handling)
scipy>=1.7 (interpolation, optimization)
numpy>=1.20 (array operations)
matplotlib>=3.3 (diagnostic plots)
Quick Start
python
Copy code

from mews import measure_equivalent_width

# Measure Hδ absorption in MaNGA datacube
result = measure_equivalent_width(
    file_path='manga-7443-12701.fits',
    file_type='fits',
    line_name='Hdelta',
    continuum_method='improved'  # 5th-degree polynomial
)

print(f"Equivalent Width: {result['ew']:.2f} ± {result['ew_err']:.2f} Å")

Usage Examples
Single Galaxy Analysis
python
Copy code

# Measure multiple lines in one spectrum
lines = ['Hdelta', 'Hgamma', 'Hbeta', 'Halpha', 'OII3727']

for line in lines:
    result = measure_equivalent_width(
        file_path='spectrum.json',
        file_type='json',
        line_name=line,
        continuum_method='improved'
    )
    print(f"{line}: {result['ew']:.2f} Å")
Batch Processing
python
Copy code

# Process entire catalog
import pandas as pd

catalog = pd.read_csv('e_a_catalog.csv')
results = []

for idx, galaxy in catalog.iterrows():
    ew_result = measure_equivalent_width(
        file_path=f"data/{galaxy['plateifu']}.fits",
        file_type='fits',
        line_name='Hdelta'
    )
    results.append(ew_result)

catalog['Hdelta_EW'] = [r['ew'] for r in results]
Continuum Fitting Methods
Improved (Recommended)
5th-degree polynomial with wavelength normalization and iterative sigma clipping (2.5σ). Best for most galaxy spectra. This method was empirically validated to reproduce manual PyRAF measurements.

Segmented
Divides spectrum into 8-10 overlapping segments, fits 2nd-degree polynomials independently. Good for complex continuum shapes.

Robust
Iterative outlier rejection using MAD (Median Absolute Deviation) statistics. Stable for noisy spectra.

Median Filter
Combines median filtering with polynomial fitting. Effective for high-noise data.

Validation & Diagnostics
MEWS produces comprehensive diagnostic plots for every measurement:

Full spectrum with fitted continuum
Local continuum fit around target line
Line profile fit (Gaussian or Voigt)
Residuals analysis
Component breakdown
Out-of-range/failed fit warnings
Example output: Visual verification at each step ensures measurement quality comparable to manual inspection.

Scientific Applications
Post-Starburst Galaxy Classification: E+A identification via Balmer absorption (Greene et al. 2021, 2026)
Stellar Population Studies: Age-sensitive absorption features
Emission Line Analysis: Star formation tracers (sign convention: absorption = positive EW)
AGN Diagnostics: [O III], [N II], [S II] measurements for BPT diagrams
Related Projects
E+A Galaxy Catalog: 183 spatially-resolved post-starburst galaxies measured with MEWS
MOONJAM: Modular analysis suite for MaNGA IFS data (in development)
Citation
If you use MEWS in your research, please cite:

bibtex
Copy code

@article{Greene2026,
  author = {Greene, Olivia A. and others},
  title = {A Complete Catalog of Post-starburst, E+A Galaxies in SDSS-IV MaNGA (MPL-11): 
           A Citizen Science Approach to Spectrophotometric Classification \& 
           the Automation of Equivalent Width Measurements},
  journal = {The Astrophysical Journal},
  year = {2026},
  note = {In Preparation}
}

@article{Greene2021,
  author = {Greene, Olivia A. and others},
  title = {Refining the E+A Galaxy: A Spatially Resolved Spectrophotometric Sample 
           of Nearby Post-starburst Systems in SDSS-IV MaNGA (MPL-5)},
  journal = {The Astrophysical Journal},
  volume = {910},
  pages = {162},
  year = {2021},
  doi = {10.3847/1538-4357/abe4d4}
}
Development History
MEWS evolved through three phases:

Manual PyRAF measurements (MPL-5 sample, 30 galaxies) - Identified systematic errors in SDSS pipeline
Jupyter Notebook prototype (interactive verification) - Developed adaptive algorithms
Full Python pipeline (MPL-11 scale, 579+ galaxies) - Production-ready automation
Development driven by need to maintain visual inspection quality at scale while addressing systematic measurement biases in automated catalogs.

Technical Details
Sign Convention
Following SDSS standards: absorption features have positive EW; emission features have negative EW

Wavelength Coverage
Optimized for optical spectroscopy (3500-10000 Å), with particular focus on:

Hydrogen Balmer series (Hα, Hβ, Hγ, Hδ)
Forbidden emission lines ([O II], [O III], [N II], [S II])
4000 Å break region
Error Estimation
Errors calculated from:

Continuum fit uncertainties
Line profile fit residuals
Local noise estimates
Contributing
Contributions welcome! Areas for development:

Additional continuum fitting algorithms
Machine learning-based continuum placement
Extended wavelength coverage (UV, NIR)
Automated outlier detection improvements
Unit tests and continuous integration
Please open an issue or pull request on GitHub.

Troubleshooting
Common issues:

"Line out of wavelength range": Check spectral coverage; some SDSS fibers have truncated blue/red ends
"Continuum fit failed": Try different continuum method (e.g., 'robust' for noisy data)
Unrealistic EW values: Review diagnostic plots; may indicate poor continuum placement or blended lines

Acknowledgments
This work was supported by:

Vanderbilt University Department of Physics & Astronomy
Fisk-Vanderbilt Bridge Program
SDSS-IV MaNGA Survey collaboration
27 citizen scientists who contributed to visual classification validation
Special thanks to:

Dr. Kelly Holley-Bockelmann (Advisor)
Dr. Charles T. Liu (Co-Advisor and Co-Chair)
License
MIT License

Copyright (c) 2026 Olivia A. Greene

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Contact
Olivia A. Greene, PhD
Astrophysicist | Pipeline Developer
Email: oliviaallegragreene@gmail.com
Website: https://galaxygreene.com
GitHub: @InfinitelyCurious

Developed at Vanderbilt University (2020-2026)
Part of dissertation: "Seeing What Is, What Was, What Could Be, What Must Not: Refining, Cataloging, and Investigating A Complete, Spatially Resolved Spectrophotometric Sample of Nearby Post-Starburst E+A Galaxies in SDSS-IV MaNGA"



