"""
MEASURING EQUIVALENT WIDTH in SPECTRA (MEWS)
=============================================

Automated pipeline for accurate equivalent width measurements in galaxy spectra,
with built-in validation and adaptive continuum fitting algorithms.

================================================================================
AUTHOR & CONTACT
================================================================================
Olivia A. Greene, PhD
Vanderbilt University (2020-2026)
Astrophysicist | Pipeline Developer

Email: oliviaallegragreene@gmail.com
Website: https://galaxygreene.com
GitHub: @InfinitelyCurious

================================================================================
PRIMARY CITATION (DISSERTATION)
================================================================================
Greene, O. A. (2026). "Seeing What Is, What Was, What Could Be, What Must Not: 
Refining, Cataloging, and Investigating A Complete, Spatially Resolved 
Spectrophotometric Sample of Nearby Post-Starburst E+A Galaxies in SDSS-IV MaNGA." 
PhD Dissertation, Vanderbilt University. 300+ pages.

Advisors: 
- Dr. Kelly Holley-Bockelmann (Vanderbilt University)
- Dr. Charles T. Liu (CUNY College of Staten Island / American Museum of Natural History)

This code implements the methodology from Chapter 3.3.4 Measuring Equivalent Width in Spectra (MEWS)

================================================================================
PUBLICATIONS BY AUTHOR
================================================================================

5. Greene, O. A., Dungee, R., Schonhut-Stasik, J., & Oldham, L. (2025). 
   "Neurodivergent in astronomy: the early-career researcher edition." 
   Nature Astronomy, 9, 1754-1757. [First-author interview feature]

4. Greene, O. A., et al. (2026). "A Complete Catalog of Post-starburst, 
   E+A Galaxies in SDSS-IV MaNGA (MPL-11): A Citizen Science Approach to 
   Spectrophotometric Classification & the Automation of Equivalent Width 
   Measurements." The Astrophysical Journal (In Preparation).

3. Ludwig, E., Iyer, K., Liu, C. T., Grwiser, E., & Greene, O. A. (2025). 
   "Comparing Star Formation Histories and Evolutionary Pathways of 
   Post-Starburst and E+A Galaxies in TNG50 and SDSS-IV MaNGA." 
   The Astrophysical Journal, 989(1), 87.

2. Greene, O. A., et al. (2021). "Refining the E+A Galaxy: A Spatially 
   Resolved Spectrophotometric Sample of Nearby Post-starburst Systems in 
   SDSS-IV MaNGA (MPL-5)." The Astrophysical Journal, 910, 162.
   DOI: 10.3847/1538-4357/abe4d0

1. Marinelli, M., Greene, O. A., Riffel, R. A., Rowlands, K., & Liu, C. 
   (2020). "SDSS-IV MaNGA: A Dwarf E+A Galaxy Quenched by AGN Feedback." 
   Research Notes of the AAS, 4, 110.

================================================================================
LICENSE
================================================================================
MIT License - See LICENSE file for details

Copyright (c) 2026 Olivia A. Greene

================================================================================
OVERVIEW
================================================================================

This script measures equivalent widths (EW) of spectral lines from astronomical 
spectra. It supports both JSON and FITS file formats and includes multiple 
continuum fitting methods.

Developed for post-starburst galaxy research; validated against manual PyRAF 
measurements for 183 E+A galaxies in SDSS-IV MaNGA DR17.

KEY FEATURES:
- 4 continuum fitting algorithms (improved, segmented, robust, median)
- Adaptive line fitting with Gaussian profiles
- Built-in validation and diagnostic plots
- CLI interface with preset configurations
- Batch processing for large catalogs

INSTALLATION REQUIREMENTS:
- Python 3.8+
- numpy>=1.20, pandas>=1.3, matplotlib>=3.3
- astropy>=5.0, specutils>=1.0, scipy>=1.7

See requirements.txt for complete dependencies.

USAGE EXAMPLES:

1. Process JSON files with preset configuration:
   python MEWS_Measuring_Equivalent_Width_in_Spectra.py --preset example_json --no-plots

2. Process FITS files with custom settings:
   python MEWS_Measuring_Equivalent_Width_in_Spectra.py --input-type fits 
   --folder-path /path/to/fits/files --redshift-csv /path/to/redshifts.csv 
   --output-csv results.csv

3. Use specific continuum method:
   python MEWS_Measuring_Equivalent_Width_in_Spectra.py --preset example_fits 
   --continuum-method robust

SUPPORTED SPECTRAL LINES:
- H-alpha (6564.61 Å) - absorption
- [NII] (6549.86, 6585.2 Å) - emission  
- [OIII] (5008.24 Å) - emission
- H-beta (4862.68 Å) - absorption
- H-gamma (4341.68 Å) - absorption
- H-delta (4102.89 Å) - absorption (strong E+A diagnostic)
- [OII] (3729.875 Å) - emission
- [SII] (6718.295, 6732.674 Å) - emission

CONTINUUM METHODS:
- improved: 5th-degree polynomial with outlier rejection (default, validated)
- segmented: Fit continuum in segments with spline interpolation
- robust: Iterative outlier rejection with MAD statistics
- median: Median filtering approach

OUTPUT:
- CSV file with EW measurements for each galaxy
- Optional diagnostic plots showing continuum fits and line profiles

For more examples, see the examples/ directory.

================================================================================
"""

# IMPORTS --------------------------------------------------------------------------------------------
import os
import json
import pandas as pd
import numpy as np
import warnings
import matplotlib.pyplot as plt
from glob import glob
from astropy import units as u
from astropy import constants
from astropy.io import fits
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum
from astropy.modeling.models import Linear1D, Gaussian1D
from astropy.visualization import quantity_support
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter, SimplexLSQFitter, LMLSQFitter
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from scipy.signal import medfilt
from scipy.interpolate import splrep, splev

# Default configuration
PLOTS_DIR = "EW_Measurements_Combined"  # Change this to your preferred directory name
CONTINUUM_METHOD = "improved"  # Options: "improved", "segmented", "robust", "median"

quantity_support()

# SPECTRAL LINE DEFINITIONS --------------------------------------------------------------------------------------------
# Rest wavelengths for spectral lines (in Angstroms)
H_alpha = 6564.61 * u.AA  # Hydrogen alpha (absorption)
NII_1 = 6549.86 * u.AA   # Nitrogen II doublet (emission)
NII_2 = 6585.2 * u.AA    # Nitrogen II doublet (emission)
OIII = 5008.24 * u.AA    # Oxygen III (emission)
H_beta = 4862.68 * u.AA  # Hydrogen beta (absorption)
H_gamma = 4341.68 * u.AA # Hydrogen gamma (absorption)
H_delta = 4102.89 * u.AA # Hydrogen delta (absorption)
OII = 3729.875 * u.AA    # Oxygen II (emission)
SII_1 = 6718.295 * u.AA  # Sulfur II doublet (emission)
SII_2 = 6732.674 * u.AA  # Sulfur II doublet (emission)

# Line configuration arrays
line_names = ['H_alpha', 'NII_1', 'NII_2', 'OIII', 'H_beta', 'H_gamma', 'H_delta', 'OII', 'SII_1', 'SII_2']
lines_to_measure = [H_alpha, NII_1, NII_2, OIII, H_beta, H_gamma, H_delta, OII, SII_1, SII_2]
line_values = [line.value for line in lines_to_measure]

# Initial width guesses for line fitting (in Angstroms)
# Optimized values for each specific line
width_guess = [5, 2, 5, 5, 10, 3, 3, 10, 10, 10]

# Line types: True = absorption, False = emission
line_absorption = [True, False, False, False, True, True, True, False, False, False]

# Custom flux unit for JSON files (erg/s/cm²/Å)
F_lambda_unit = u.erg / u.s / u.cm ** 2 / u.AA

print("Line definitions loaded successfully")

# PRESET CONFIGURATIONS --------------------------------------------------------------------------------------------
# Example configurations for different data types and sources
# Replace these paths with your actual file locations
PRESETS = {
    "example_json": {
        "input_type": "json",
        "folder_path": "YOUR/PATH/TO/JSON/SPECTRAL/FILES",
        "redshift_csv": "YOUR/PATH/TO/redshift_catalog.csv",
        "output_csv": "JSON_EW_measurements.csv",
        "continuum_method": "improved",
        "plots_dir": "EW_Plots_JSON"
    },
    "example_json_alt": {
        "input_type": "json",
        "folder_path": "YOUR/PATH/TO/ALTERNATIVE/JSON/FILES",
        "redshift_csv": "YOUR/PATH/TO/redshift_catalog.csv",
        "output_csv": "JSON_alt_EW_measurements.csv",
        "continuum_method": "improved",
        "plots_dir": "EW_Plots_JSON_Alt"
    },
    "example_fits": {
        "input_type": "fits",
        "folder_path": "YOUR/PATH/TO/FITS/CUBE/FOLDERS",
        "redshift_csv": "YOUR/PATH/TO/redshift_catalog.csv",
        "output_csv": "FITS_EW_measurements.csv",
        "continuum_method": "improved",
        "plots_dir": "EW_Plots_FITS"
    },
    "example_fits_alt": {
        "input_type": "fits",
        "folder_path": "YOUR/PATH/TO/FITS/CUBE/FOLDERS",
        "redshift_csv": "YOUR/PATH/TO/redshift_catalog.csv",
        "output_csv": "FITS_alt_EW_measurements.csv",
        "continuum_method": "improved",
        "plots_dir": "EW_Plots_FITS_Alt"
    }
}

# UTILITY FUNCTIONS --------------------------------------------------------------------------------------------

def validate_line_in_range(wavelength_corrected_value, line_wavelength, buffer=20):
    """
    Verify that a spectral line is within the measured wavelength range.
    
    Parameters:
    -----------
    wavelength_corrected_value : array
        Redshift-corrected wavelength array
    line_wavelength : astropy.units.Quantity
        Rest wavelength of the spectral line
    buffer : float
        Buffer zone in Angstroms from spectrum edges
        
    Returns:
    --------
    bool : True if line is within range, False otherwise
    """
    min_wave = np.min(wavelength_corrected_value) + buffer
    max_wave = np.max(wavelength_corrected_value) - buffer

    if line_wavelength.value < min_wave or line_wavelength.value > max_wave:
        print(f"WARNING: Line at {line_wavelength.value}Å is outside the spectrum range "
              f"({min_wave:.1f}Å - {max_wave:.1f}Å)")
        return False
    return True

def get_window_size(line_center):
    """
    Return appropriate spectral window size based on line wavelength.
    
    Parameters:
    -----------
    line_center : float
        Central wavelength of the line in Angstroms
        
    Returns:
    --------
    float : Window size in Angstroms
    """
    if line_center >= 6500:  # Red region (H-alpha, SII, NII)
        return 100
    elif line_center >= 4800:  # Mid region (H-beta, OIII)
        return 80
    else:  # Blue region (H-gamma, H-delta, OII)
        return 50

# CONTINUUM FITTING METHODS --------------------------------------------------------------------------------------------

def improved_continuum_fit(spectrum):
    """
    Polynomial continuum fitting with outlier rejection for galaxy spectra.
    
    This method fits a 5th-degree polynomial to the spectrum, applies sigma clipping
    to remove emission/absorption lines, then refits for a clean continuum.
    
    Parameters:
    -----------
    spectrum : Spectrum1D
        Input spectrum object
        
    Returns:
    --------
    function : Continuum function that can be evaluated at any wavelength
    """
    # Extract data
    x = spectrum.spectral_axis.value
    y = spectrum.flux.value

    # Normalize wavelengths to prevent numerical issues
    x_mean = np.mean(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_mean) / (x_range / 2)

    # Fit 5th-degree polynomial
    model = models.Polynomial1D(degree=5)
    fitter = fitting.LinearLSQFitter()
    model_fitted = fitter(model, x_norm, y)

    # Apply sigma clipping to mask emission/absorption lines
    residuals = y - model_fitted(x_norm)
    mask = ~sigma_clip(residuals, sigma=2.5).mask

    # Refit using only continuum points
    model_fitted = fitter(model, x_norm[mask], y[mask])

    # Return function that handles unit conversion
    return lambda x: model_fitted((x.value - x_mean) / (x_range / 2)) * spectrum.flux.unit

def segmented_continuum_fit(spectrum, segments=8):
    """
    Fit continuum by dividing spectrum into segments with spline interpolation.
    
    This method breaks the spectrum into overlapping segments, fits low-degree
    polynomials to each segment with sigma clipping, then creates a smooth
    spline through all segment fits.
    
    Parameters:
    -----------
    spectrum : Spectrum1D
        Input spectrum
    segments : int
        Number of segments to divide spectrum into
        
    Returns:
    --------
    function : Spline-interpolated continuum function
    """
    x = spectrum.spectral_axis.value
    y = spectrum.flux.value

    # Divide wavelength range into overlapping segments
    x_min, x_max = np.min(x), np.max(x)
    segment_size = (x_max - x_min) / segments

    x_full = np.array([])
    y_cont_full = np.array([])

    for i in range(segments):
        # Define overlapping segment bounds (15% overlap)
        seg_min = x_min + i * segment_size - segment_size * 0.15
        seg_max = x_min + (i + 1) * segment_size + segment_size * 0.15

        # Filter data for this segment
        mask = (x >= seg_min) & (x <= seg_max)
        if not np.any(mask) or np.sum(mask) < 10:
            continue

        x_seg = x[mask]
        y_seg = y[mask]

        # Normalize for numerical stability
        x_mean = np.mean(x_seg)
        x_range = np.max(x_seg) - np.min(x_seg)
        x_norm = (x_seg - x_mean) / (x_range / 2)

        # Fit low-degree polynomial with aggressive sigma clipping
        model = models.Polynomial1D(degree=2)
        fitter = fitting.LinearLSQFitter()
        model_fitted = fitter(model, x_norm, y_seg)

        # Remove emission/absorption lines and refit
        residuals = y_seg - model_fitted(x_norm)
        clean_mask = ~sigma_clip(residuals, sigma=1.5).mask

        if np.sum(clean_mask) > 5:
            model_fitted = fitter(model, x_norm[clean_mask], y_seg[clean_mask])

        # Store segment results
        y_cont_seg = model_fitted(x_norm)
        x_full = np.append(x_full, x_seg)
        y_cont_full = np.append(y_cont_full, y_cont_seg)

    # Sort by wavelength and create smoothing spline
    sort_idx = np.argsort(x_full)
    x_full = x_full[sort_idx]
    y_cont_full = y_cont_full[sort_idx]

    spline = splrep(x_full, y_cont_full, s=0.5*len(x_full))

    return lambda wavelength: splev(wavelength.value, spline) * spectrum.flux.unit

def robust_continuum_fit(spectrum):
    """
    Robust continuum fitting with iterative outlier rejection using MAD statistics.
    
    This method uses median absolute deviation (MAD) for robust outlier detection
    and iteratively refits the continuum while rejecting spectral features.
    
    Parameters:
    -----------
    spectrum : Spectrum1D
        Input spectrum
        
    Returns:
    --------
    function : Continuum function with robust outlier rejection
    """
    x = spectrum.spectral_axis.value
    y = spectrum.flux.value

    # Normalize wavelengths for numerical stability
    x_mean = np.mean(x)
    x_range = np.max(x) - np.min(x)
    x_norm = (x - x_mean) / (x_range / 2)

    # Initial medium-degree polynomial fit
    model = models.Polynomial1D(degree=4)
    fitter = fitting.LinearLSQFitter()
    initial_fit = fitter(model, x_norm, y)

    # Iterative outlier rejection using MAD statistics (3 iterations)
    model_fitted = initial_fit
    residuals = y - initial_fit(x_norm)

    for iteration in range(3):
        # Calculate robust statistics using MAD
        med_resid = np.median(residuals)
        mad = np.median(np.abs(residuals - med_resid))

        # Identify inlier points (within 2.5 MAD from median)
        mask = np.abs(residuals - med_resid) < 2.5 * mad

        # Ensure we don't lose too many points
        if np.sum(mask) < len(x) * 0.5:
            break

        # Refit with inliers only
        model_fitted = fitter(model, x_norm[mask], y[mask])
        residuals = y - model_fitted(x_norm)

    return lambda wavelength: model_fitted((wavelength.value - x_mean) / (x_range / 2)) * spectrum.flux.unit

def median_filter_continuum(spectrum, window_size=101, poly_degree=3):
    """
    Median filtering approach for continuum estimation.
    
    Applies median filtering to smooth out spectral features, then fits
    a polynomial to the filtered result for a smooth continuum.
    
    Parameters:
    -----------
    spectrum : Spectrum1D
        Input spectrum
    window_size : int
        Size of median filter window (must be odd)
    poly_degree : int
        Degree of polynomial fit to filtered data
        
    Returns:
    --------
    function : Median-filtered continuum function
    """
    x = spectrum.spectral_axis.value
    y = spectrum.flux.value

    # Apply median filter to suppress spectral features
    y_med = medfilt(y, window_size)

    # Fit polynomial to smoothed data for final continuum
    coeff = np.polyfit(x, y_med, poly_degree)

    return lambda wavelength: np.polyval(coeff, wavelength.value) * spectrum.flux.unit

def get_continuum_fit(spectrum, method=CONTINUUM_METHOD, line_center=None, width_guess=5.0):
    """
    Select and apply appropriate continuum fitting method with optional line masking.
    
    Parameters:
    -----------
    spectrum : Spectrum1D
        Input spectrum
    method : str
        Continuum fitting method to use
    line_center : float, optional
        Central wavelength for line masking during continuum fit
    width_guess : float, optional
        Line width estimate for masking
        
    Returns:
    --------
    function : Selected continuum fitting function
    """
    # Apply line masking if line center is provided
    if line_center is not None:
        x = spectrum.spectral_axis.value
        y = spectrum.flux.value

        # Mask line region (exclude points within 3x the width_guess)
        line_mask = np.abs(x - line_center) > width_guess * 3.0

        # Only apply mask if sufficient points remain for continuum fitting
        if np.sum(line_mask) >= 10:
            masked_spectrum = Spectrum1D(
                flux=spectrum.flux[line_mask],
                spectral_axis=spectrum.spectral_axis[line_mask]
            )
            spectrum = masked_spectrum

    # Select continuum fitting method
    if method == "improved":
        return improved_continuum_fit(spectrum)
    elif method == "segmented":
        return segmented_continuum_fit(spectrum)
    elif method == "robust":
        return robust_continuum_fit(spectrum)
    elif method == "median":
        return median_filter_continuum(spectrum)
    else:
        # Default to improved method
        return improved_continuum_fit(spectrum)

# DATA I/O FUNCTIONS --------------------------------------------------------------------------------------------

def read_trace_keys_from_json(json_path):
    """
    Read available trace keys from a JSON spectral file.
    
    Parameters:
    -----------
    json_path : str
        Path to JSON file containing spectral traces
        
    Returns:
    --------
    list : Available trace keys in the JSON file
    """
    with open(json_path, 'r') as file:
        data = json.load(file)
    return list(data['traces'].keys())

def read_spectrum_from_json(json_path, trace_key):
    """
    Extract wavelength and flux data from JSON file for a specific trace.
    
    Parameters:
    -----------
    json_path : str
        Path to JSON spectral file
    trace_key : str
        Specific trace identifier to extract
        
    Returns:
    --------
    tuple : (wavelength, flux) arrays with appropriate units
    """
    with open(json_path, 'r') as file:
        data = json.load(file)

    trace_data = data['traces'][trace_key]

    # Process wavelength and flux with proper units and scaling
    wavelength = np.array(trace_data['wavelength']) * u.Angstrom
    flux = np.array(trace_data['flux']) * F_lambda_unit * 10**-17  # Apply scaling factor
    return wavelength, flux

# LINE FITTING AND EW MEASUREMENT --------------------------------------------------------------------------------------------

def fit_line_profile(filtered_wavelength_value, spec_normalized_value, line_wavelength, width_guess, absorption=True):
    """
    Fit spectral line profile with automatic Gaussian/Voigt model selection.
    
    This function automatically selects between Gaussian and Voigt profiles based
    on line depth/strength. Deep lines use Voigt profiles for better modeling of
    line wings, while normal lines use simpler Gaussian profiles.
    
    Parameters:
    -----------
    filtered_wavelength_value : array
        Wavelength array for the line region
    spec_normalized_value : array  
        Continuum-normalized flux array
    line_wavelength : astropy.units.Quantity
        Rest wavelength of the line
    width_guess : float
        Initial width estimate in Angstroms
    absorption : bool
        True for absorption lines, False for emission
        
    Returns:
    --------
    tuple : (fitted_model, is_deep_line) where fitted_model is the best-fit model
            and is_deep_line indicates if Voigt profile was used
    """
    from astropy.modeling.models import Voigt1D

    # Identify line region around expected position
    line_region = np.abs(filtered_wavelength_value - line_wavelength.value) < width_guess * 2

    if not np.any(line_region):
        raise ValueError("No data points in line region")

    # Estimate line strength and initial amplitude
    if absorption:
        line_value = np.min(spec_normalized_value[line_region])
        amp_guess = -(1.0 - line_value) * 1.5  # Negative for absorption
    else:
        line_value = np.max(spec_normalized_value[line_region])
        amp_guess = (line_value - 1.0) * 1.5  # Positive for emission

    # Determine if this is a deep line requiring Voigt profile
    is_deep_line = (absorption and line_value < 0.8) or (not absorption and line_value > 1.2)

    # Set up appropriate model based on line depth
    if is_deep_line:
        # Use Voigt for deep lines with broader wings
        fwhm_g = width_guess * 2.355  # Gaussian FWHM
        fwhm_l = width_guess * 1.5    # Lorentzian FWHM

        # Increase amplitude for very deep absorption lines
        amp_factor = 2.0 if (absorption and line_value < 0.7) else 1.5

        line_model = models.Linear1D(slope=0, intercept=1.0) + \
                    Voigt1D(x_0=line_wavelength.value,
                           amplitude_L=amp_guess * amp_factor,
                           fwhm_L=fwhm_l, fwhm_G=fwhm_g)
    else:
        # Standard Gaussian for normal lines
        line_model = models.Linear1D(slope=0, intercept=1.0) + \
                   models.Gaussian1D(amplitude=amp_guess,
                                    mean=line_wavelength.value,
                                    stddev=width_guess)

    # Fix continuum parameters and line center
    line_model[0].intercept.fixed = True
    line_model[0].slope.fixed = True

    # Fix line center (wavelength position)
    if hasattr(line_model[1], 'mean'):
        line_model[1].mean.fixed = True
    elif hasattr(line_model[1], 'x_0'):
        line_model[1].x_0.fixed = True

    # Try multiple fitters with different algorithms
    fitters = [
        (fitting.LevMarLSQFitter(), {'maxiter': 2000}),
        (fitting.SimplexLSQFitter(), {}),
        (fitting.LMLSQFitter(), {})
    ]

    # Attempt fitting with each fitter until one succeeds
    for fitter, kwargs in fitters:
        try:
            fitted_model = fitter(line_model, filtered_wavelength_value, spec_normalized_value, **kwargs)

            # Quality check - ensure reasonable amplitude
            if hasattr(fitted_model[1], 'amplitude'):
                if abs(fitted_model[1].amplitude.value) < 0.001:
                    continue  # Amplitude too small, try next fitter

            return fitted_model, is_deep_line

        except Exception:
            continue  # Fitting failed, try next fitter

    # Create fallback model if all fitters fail
    amp = min(-(1.0 - line_value) * 2, -0.9) if absorption else max((line_value - 1.0) * 2, 0.9)

    if is_deep_line:
        basic_model = models.Linear1D(slope=0, intercept=1.0) + \
                     Voigt1D(x_0=line_wavelength.value,
                            amplitude_L=amp,
                            fwhm_L=width_guess*3,
                            fwhm_G=width_guess*2.355)
    else:
        basic_model = models.Linear1D(slope=0, intercept=1.0) + \
                     models.Gaussian1D(amplitude=amp,
                                     mean=line_wavelength.value,
                                     stddev=width_guess)

    print(f"Warning: Used fallback model for {line_wavelength.value}Å")
    return basic_model, is_deep_line

def measure_ew(wavelength, flux, z, line_wavelength, width_guess, absorption=True, plot=False, save_dir=PLOTS_DIR, plot_full=False):
    """
    Measure the Equivalent Width (EW) for a spectral line.
    
    This is the main function that performs the complete EW measurement process:
    1. Applies redshift correction to observed wavelengths
    2. Extracts the spectral region around the line
    3. Fits the local continuum 
    4. Normalizes the spectrum by the continuum
    5. Fits the line profile (Gaussian or Voigt)
    6. Calculates the equivalent width
    7. Creates diagnostic plots (optional)
    
    Parameters:
    -----------
    wavelength : astropy.units.Quantity
        Observed wavelength array
    flux : astropy.units.Quantity  
        Observed flux array
    z : float
        Redshift of the galaxy
    line_wavelength : astropy.units.Quantity
        Rest wavelength of the line to measure
    width_guess : float
        Initial estimate of line width in Angstroms
    absorption : bool
        True for absorption lines, False for emission lines
    plot : bool
        Whether to create diagnostic plots
    save_dir : str
        Directory to save plots
    plot_full : bool
        Whether to plot the full spectrum continuum fit
        
    Returns:
    --------
    astropy.units.Quantity : Equivalent width in Angstroms
                            Returns NaN if measurement fails
    """
    os.makedirs(save_dir, exist_ok=True)

    # Get appropriate spectral window size based on wavelength
    window_size = get_window_size(line_wavelength.value)
    print(f"Using window size {window_size}Å for line at {line_wavelength.value}Å")

    # Apply redshift correction to move to rest frame
    wavelength_value = wavelength.value
    wavelength_corrected_value = wavelength_value / (1 + z)
    wavelength_corrected = wavelength_corrected_value * u.AA

    flux_value = flux.value

    # Validate that line is within the spectrum's wavelength coverage
    if not validate_line_in_range(wavelength_corrected_value, line_wavelength):
        print(f"Skipping measurement for {line_wavelength.value}Å (outside spectrum range)")
        if plot:
            # Create diagnostic plot for out-of-range lines
            plt.figure(figsize=(10, 6))
            plt.plot(wavelength_corrected_value, flux_value, label='Spectrum')
            plt.axvline(x=line_wavelength.value, color='red', linestyle='--',
                       label=f'Line Position ({line_wavelength.value:.1f}Å) - OUTSIDE RANGE')
            plt.text(0.05, 0.95,
                    f"WARNING: Line at {line_wavelength.value:.1f}Å\n"
                    f"is outside spectrum range\n"
                    f"({np.min(wavelength_corrected_value):.1f}Å - {np.max(wavelength_corrected_value):.1f}Å)",
                    transform=plt.gca().transAxes,
                    fontsize=10, bbox=dict(facecolor='red', alpha=0.2))
            plt.xlabel('Wavelength (Å)')
            plt.ylabel('Flux')
            plt.title(f"Line Outside Range: {line_wavelength.value}Å")
            plt.legend()

            filename = f"{save_dir}/outside_range_{int(line_wavelength.value)}_angstrom.png"
            plt.savefig(filename)
            plt.close()

        return np.nan * u.AA

    # Plot full spectrum continuum fit (for first line of each galaxy only)
    if plot and plot_full and line_wavelength == lines_to_measure[0]:
        full_spectrum_flux = flux_value * flux.unit
        full_spectrum_wave = wavelength_corrected_value * u.AA
        full_spectrum = Spectrum1D(flux=full_spectrum_flux, spectral_axis=full_spectrum_wave)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                # Use segmented method with more segments for full spectrum
                full_cont_fit = segmented_continuum_fit(full_spectrum, segments=10)
                full_continuum_with_units = full_cont_fit(full_spectrum_wave)
                full_continuum_value = full_continuum_with_units.value

                # Create full spectrum diagnostic plot
                plt.figure(figsize=(12, 6))
                plt.plot(wavelength_corrected_value, flux_value,
                        label="Spectrum", alpha=0.7)
                plt.plot(wavelength_corrected_value, full_continuum_value,
                        label="Segmented Continuum Fit", color='red', linewidth=2)

                # Mark positions of all measured spectral lines
                for line_name, line_wave, is_abs in zip(line_names, lines_to_measure, line_absorption):
                    if validate_line_in_range(wavelength_corrected_value, line_wave):
                        plt.axvline(x=line_wave.value, color='green', linestyle='--', alpha=0.7)
                        y_pos = plt.ylim()[1] * 0.9 if plt.ylim()[1] > 0 else plt.ylim()[1] * 0.1
                        plt.text(line_wave.value, y_pos, f"{line_name}",
                                rotation=90, verticalalignment='top')

                plt.xlabel('Wavelength (Å)')
                plt.ylabel('Flux')
                plt.title('Full Spectrum with Segmented Continuum Fit')
                plt.legend()

                filename = f"{save_dir}/full_spectrum_continuum.png"
                plt.savefig(filename)
                plt.close()
                print(f"Saved full spectrum plot to {filename}")
            except Exception as e:
                print(f"Error creating full spectrum plot: {e}")

    # Extract spectral region around the line
    line_filter = np.abs(wavelength_corrected_value - line_wavelength.value) < window_size
    filtered_wavelength_value = wavelength_corrected_value[line_filter]
    filtered_flux_value = flux_value[line_filter]

    # Ensure adequate data points for analysis
    if len(filtered_wavelength_value) < 5:
        print(f"Warning: Too few points in window for {line_wavelength.value}Å. Increasing window size.")
        window_size *= 1.5
        line_filter = np.abs(wavelength_corrected_value - line_wavelength.value) < window_size
        filtered_wavelength_value = wavelength_corrected_value[line_filter]
        filtered_flux_value = flux_value[line_filter]

    # Create spectrum object with proper units
    filtered_wavelength = filtered_wavelength_value * u.AA
    filtered_flux = filtered_flux_value * flux.unit
    spectrum = Spectrum1D(flux=filtered_flux, spectral_axis=filtered_wavelength)

    # Fit local continuum around the line
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            # Pass line center for improved continuum fitting with line masking
            g1_fit = get_continuum_fit(spectrum, method=CONTINUUM_METHOD,
                                     line_center=line_wavelength.value,
                                     width_guess=width_guess)
            y_continuum_fitted_with_units = g1_fit(filtered_wavelength)
            y_continuum_fitted_value = y_continuum_fitted_with_units.value
        except Exception as e:
            print(f"Continuum fitting failed: {e}. Using median flux.")
            # Fallback to simple median continuum
            y_continuum_fitted_value = np.ones_like(filtered_wavelength_value) * np.median(filtered_flux_value)

    # Create continuum diagnostic plot
    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(filtered_wavelength_value, filtered_flux_value, label='Spectrum')
        plt.plot(filtered_wavelength_value, y_continuum_fitted_value, color='purple', label='Continuum Fit')
        plt.axvline(x=line_wavelength.value, color='green', linestyle='--',
                   label=f'Line Position ({line_wavelength.value:.1f}Å)')
        plt.text(0.05, 0.95, f"Window size: {window_size}Å", transform=plt.gca().transAxes,
                 fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux')
        plt.legend()
        plt.title(f"Continuum Fit for {line_wavelength.value} Angstrom")

        filename = f"{save_dir}/continuum_fit_{int(line_wavelength.value)}_angstrom.png"
        plt.savefig(filename)
        plt.close()
        print(f"Saved continuum plot to {filename}")

    # Normalize spectrum by dividing by continuum
    spec_normalized_value = filtered_flux_value / y_continuum_fitted_value

    # Fit line profile and calculate equivalent width
    try:
        # Fit line profile using enhanced function
        fitted_model, is_deep_line = fit_line_profile(
            filtered_wavelength_value,
            spec_normalized_value,
            line_wavelength,
            width_guess,
            absorption=absorption
        )

        # Extract model components
        continuum_model, profile_model = fitted_model
        profile_values = profile_model(filtered_wavelength_value)
        model_values = fitted_model(filtered_wavelength_value)

        # Calculate equivalent width from line profile integral
        delta_wavelength = np.median(np.diff(filtered_wavelength_value))
        ew_value = -np.sum(profile_values) * delta_wavelength

        # Apply correct sign convention
        # Absorption lines: EW > 0, Emission lines: EW < 0  
        if absorption and ew_value < 0:
            ew_value = -ew_value
        elif not absorption and ew_value > 0:
            ew_value = -ew_value

        ew = ew_value * u.AA

        # Get line width information for display
        if is_deep_line:
            width_info = f"FWHM(G)={profile_model.fwhm_G.value:.2f}Å, FWHM(L)={profile_model.fwhm_L.value:.2f}Å"
            model_type = "Voigt"
        else:
            width_info = f"σ={profile_model.stddev.value:.2f}Å, FWHM={profile_model.stddev.value*2.355:.2f}Å"
            model_type = "Gaussian"

        # Create diagnostic plots
        if plot:
            # Line fit plot
            plt.figure(figsize=(10, 6))
            plt.plot(filtered_wavelength_value, spec_normalized_value, label='Normalized Spectrum')
            plt.plot(filtered_wavelength_value, model_values, label='Fitted Model', color='red')
            plt.axvline(x=line_wavelength.value, color='green', linestyle='--',
                       label=f'Line Position ({line_wavelength.value:.1f}Å)')
            plt.text(0.05, 0.05,
                    f"Using {model_type} profile\n"
                    f"EW = {ew.value:.2f} Å\n"
                    f"Window size: {window_size}Å\n"
                    f"{width_info}",
                    transform=plt.gca().transAxes,
                    fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))

            plt.xlabel('Wavelength (Å)')
            plt.ylabel('Normalized Flux')
            plt.axhline(y=1.0, color='grey', linestyle='--', alpha=0.5)
            plt.legend()
            plt.title(f"Line Fit for {line_wavelength.value} Angstrom")

            filename = f"{save_dir}/line_fit_{int(line_wavelength.value)}_angstrom.png"
            plt.savefig(filename)
            plt.close()
            print(f"Saved line fit plot to {filename}")

            # Component breakdown plot
            plt.figure(figsize=(10, 6))
            plt.plot(filtered_wavelength_value, spec_normalized_value, label='Normalized Spectrum')
            plt.plot(filtered_wavelength_value, profile_values + 1.0,
                    label=f'{model_type} Component', color='green')
            plt.axhline(y=1.0, color='blue', linestyle='--', label='Continuum Level')
            plt.axvline(x=line_wavelength.value, color='green', linestyle='-.',
                       label=f'Line Position ({line_wavelength.value:.1f}Å)')
            plt.text(0.05, 0.05,
                    f"EW = {ew.value:.2f} Å\n"
                    f"Window size: {window_size}Å\n"
                    f"{width_info}",
                    transform=plt.gca().transAxes,
                    fontsize=10, bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))

            plt.xlabel('Wavelength (Å)')
            plt.ylabel('Normalized Flux')
            plt.ylim(min(0.0, min(spec_normalized_value) - 0.1), max(1.1, max(spec_normalized_value) + 0.1))
            plt.legend()
            plt.title(f"Components for {line_wavelength.value} Angstrom")

            filename = f"{save_dir}/components_{int(line_wavelength.value)}_angstrom.png"
            plt.savefig(filename)
            plt.close()
            print(f"Saved components plot to {filename}")

            # Residuals plot
            plt.figure(figsize=(10, 6))
            plt.plot(filtered_wavelength_value, spec_normalized_value, label='Normalized Spectrum', color='blue')
            residuals = spec_normalized_value - model_values
            plt.plot(filtered_wavelength_value, residuals + 1.0,
                     label='Residuals', color='orange')
            plt.axhline(y=1.0, color='grey', linestyle='--', alpha=0.5, label='Continuum/Zero Level')
            plt.xlabel('Wavelength (Å)')
            plt.ylabel('Normalized Flux')
            plt.ylim(min(0.5, min(residuals + 1.0)-0.1), max(1.5, max(spec_normalized_value)+0.1))
            plt.title(f"Residuals for {line_wavelength.value} Angstrom")
            plt.legend()

            filename = f"{save_dir}/residuals_{int(line_wavelength.value)}_angstrom.png"
            plt.savefig(filename)
            plt.close()
            print(f"Saved residuals plot to {filename}")

        print(f"{line_wavelength}: EW = {ew}")
        return ew

    except Exception as e:
        print(f"Error fitting line at {line_wavelength}: {e}")
        if plot:
            # Create diagnostic plot for failed fits
            plt.figure(figsize=(10, 6))
            plt.plot(filtered_wavelength_value, spec_normalized_value, label='Normalized Spectrum')
            plt.axhline(y=1.0, color='grey', linestyle='--', alpha=0.5)
            plt.axvline(x=line_wavelength.value, color='green', linestyle='--',
                       label=f'Line Position ({line_wavelength.value:.1f}Å)')
            plt.text(0.05, 0.95,
                    f"Window size: {window_size}Å\n"
                    f"WARNING: Fit Failed\n"
                    f"Error: {str(e)[:50]}",
                    transform=plt.gca().transAxes,
                    fontsize=10, bbox=dict(facecolor='red', alpha=0.2))
            plt.xlabel('Wavelength (Å)')
            plt.ylabel('Normalized Flux')
            plt.title(f"Failed Fit for {line_wavelength.value} Angstrom")
            plt.legend()

            filename = f"{save_dir}/failed_fit_{int(line_wavelength.value)}_angstrom.png"
            plt.savefig(filename)
            plt.close()
            print(f"Saved failed fit diagnostic plot to {filename}")

        return np.nan * u.AA

# FILE PROCESSING FUNCTIONS --------------------------------------------------------------------------------------------

def process_json_file(json_file, df, results_df, disable_plots=False):
    """
    Process a single JSON spectral file to extract EW measurements.
    
    Parameters:
    -----------
    json_file : str
        Path to JSON spectral file
    df : pandas.DataFrame  
        DataFrame containing redshift information
    results_df : pandas.DataFrame
        Existing results DataFrame to append to
    disable_plots : bool
        Whether to disable diagnostic plotting
        
    Returns:
    --------
    pandas.DataFrame : Updated results DataFrame
    """
    plate_ifu = os.path.basename(json_file).split('.')[0]

    try:
        # Look up redshift for this galaxy
        if not df.loc[df['PLATE-IFU'] == plate_ifu, 'REDSHIFT'].empty:
            z = df.loc[df['PLATE-IFU'] == plate_ifu, 'REDSHIFT'].iloc[0]
            print(f"Processing PLATE-IFU: {plate_ifu} with redshift {z}")
        else:
            print(f"No redshift found for PLATE-IFU: {plate_ifu}")
            return results_df

        # Read spectrum from JSON file
        trace_keys = read_trace_keys_from_json(json_file)
        trace_key = trace_keys[0]
        wavelength, flux = read_spectrum_from_json(json_file, trace_key)

        # Create output directory for this galaxy
        galaxy_plots_dir = f"{PLOTS_DIR}/{plate_ifu}"
        if not disable_plots:
            os.makedirs(galaxy_plots_dir, exist_ok=True)

        # Measure EW for each spectral line
        all_eqw = []
        for i, (line_name, line_wavelength, guess, is_absorption) in enumerate(zip(line_names, lines_to_measure, width_guess, line_absorption)):
            try:
                # Plot full spectrum for first line only
                plot_full_spectrum = (i == 0) and not disable_plots

                ew = measure_ew(wavelength, flux, z, line_wavelength, width_guess=guess,
                              absorption=is_absorption,
                              plot=not disable_plots,
                              plot_full=plot_full_spectrum,
                              save_dir=galaxy_plots_dir)

                if hasattr(ew, 'value'):
                    all_eqw.append(ew.value)
                else:
                    all_eqw.append(None)
            except Exception as e:
                print(f"Failed to measure EW for {line_name}: {e}")
                all_eqw.append(None)

        # Store results
        result_data = {'PLATE-IFU': plate_ifu, 'z': z}
        result_data.update({line_name: eqw for line_name, eqw in zip(line_names, all_eqw)})

        result_row_df = pd.DataFrame([result_data])
        return pd.concat([results_df, result_row_df], ignore_index=True)

    except Exception as e:
        print(f"Error processing {plate_ifu}: {e}")
        return results_df

def process_fits_file(fits_file, plate_ifu, z, results_df, disable_plots=False):
    """
    Process a single FITS cube file to extract EW measurements from central spaxel.
    
    Parameters:
    -----------
    fits_file : str
        Path to FITS cube file
    plate_ifu : str
        Galaxy identifier
    z : float
        Redshift of the galaxy
    results_df : pandas.DataFrame
        Existing results DataFrame to append to
    disable_plots : bool
        Whether to disable diagnostic plotting
        
    Returns:
    --------
    pandas.DataFrame : Updated results DataFrame
    """
    try:
        print(f"Processing PLATE-IFU: {plate_ifu} with redshift {z}")

        # Open FITS file and extract central spaxel spectrum
        HDU = fits.open(fits_file)
        flux = HDU[1].data
        wavedata_raw = HDU[6].data * u.AA

        # Get central spaxel coordinates from header
        x = int(HDU[1].header['CRPIX1'])
        y = int(HDU[1].header['CRPIX2'])

        # Extract spectrum at central spaxel
        middleflux = flux[:, x, y]

        # Create output directory
        galaxy_plots_dir = f"{PLOTS_DIR}/{plate_ifu}"
        if not disable_plots:
            os.makedirs(galaxy_plots_dir, exist_ok=True)

        # Measure EW for each spectral line
        all_eqw = []
        for i, (line_name, line_wavelength, guess, is_absorption) in enumerate(zip(line_names, lines_to_measure, width_guess, line_absorption)):
            try:
                plot_full_spectrum = (i == 0) and not disable_plots

                ew = measure_ew(wavedata_raw, middleflux * u.dimensionless_unscaled,
                              z, line_wavelength, width_guess=guess,
                              absorption=is_absorption,
                              plot=not disable_plots,
                              plot_full=plot_full_spectrum,
                              save_dir=galaxy_plots_dir)

                if hasattr(ew, 'value'):
                    all_eqw.append(ew.value)
                else:
                    all_eqw.append(np.nan)
            except Exception as e:
                print(f"Failed to measure EW for {line_name}: {e}")
                all_eqw.append(np.nan)

        HDU.close()

        # Store results with explicit column mapping
        result_data = {'PLATE-IFU': plate_ifu, 'z': z}
        for i, name in enumerate(line_names):
            if i < len(all_eqw):
                result_data[name] = all_eqw[i]
            else:
                result_data[name] = np.nan

        # Add to results DataFrame
        new_row = pd.DataFrame([result_data])
        results_df = pd.concat([results_df, new_row], ignore_index=True)

        return results_df

    except Exception as e:
        print(f"Error processing {plate_ifu}: {e}")
        return results_df

def process_json_files(folder_path, redshift_csv_path, output_csv_path, disable_plots=False):
    """
    Process all JSON files in a directory structure.
    
    Parameters:
    -----------
    folder_path : str
        Path to directory containing JSON spectral files
    redshift_csv_path : str  
        Path to CSV file with galaxy redshifts
    output_csv_path : str
        Path to save results CSV file
    disable_plots : bool
        Whether to disable diagnostic plotting
        
    Returns:
    --------
    pandas.DataFrame : Complete results DataFrame
    """
    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Find all JSON files in the directory
    json_files = glob(f'{folder_path}/*.json')
    print(f"Found {len(json_files)} JSON files in {folder_path}")

    # Load redshift catalog
    df = pd.read_csv(redshift_csv_path, index_col=False)
    results_df = pd.DataFrame(columns=['PLATE-IFU', 'z'] + line_names)

    # Process each JSON file
    for json_file in json_files:
        results_df = process_json_file(json_file, df, results_df, disable_plots)

        # Save intermediate results to prevent data loss
        results_df.to_csv(output_csv_path, index=False)
        print(f"Intermediate results saved to {output_csv_path}")

    # Final save
    results_df.to_csv(output_csv_path, index=False)
    print(f"Processing complete. Results saved to {output_csv_path}")
    return results_df

def process_fits_files(folder_path, redshift_csv_path, output_csv_path, disable_plots=False):
    """
    Process all FITS cube files in a directory structure.
    
    This function searches for FITS files in subdirectories, matches them with
    redshift information, and processes each galaxy's central spaxel spectrum.
    
    Parameters:
    -----------
    folder_path : str
        Path to directory containing FITS cube subdirectories
    redshift_csv_path : str
        Path to CSV file with galaxy redshifts  
    output_csv_path : str
        Path to save results CSV file
    disable_plots : bool
        Whether to disable diagnostic plotting
        
    Returns:
    --------
    pandas.DataFrame : Complete results DataFrame
    """
    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Search for FITS files in subdirectory structure
    filelist = glob(f'{folder_path}/*/*.fits.gz')
    print(f"Found {len(filelist)} FITS files in {folder_path}")

    if len(filelist) == 0:
        print("No FITS files found. Please check the folder path structure.")
        return pd.DataFrame()

    # Load redshift catalog
    redshift_df = pd.read_csv(redshift_csv_path, index_col=False)

    # Initialize empty results DataFrame
    results_df = pd.DataFrame(columns=['PLATE-IFU', 'z'] + line_names)
    measured_count = 0

    # Process each FITS file found
    for fits_file in filelist:
        filename = os.path.basename(fits_file)

        # Validate filename format (expected: manga-PLATE-IFU-LOGCUBE.fits.gz)
        if not filename.startswith('manga-') or not '-LOGCUBE' in filename:
            print(f"Skipping non-standard filename: {filename}")
            continue

        # Extract PLATE-IFU identifier from filename
        plate_ifu = filename.split('manga-')[1].split('-LOGCUBE')[0]

        # Look up redshift for this galaxy
        matching_rows = redshift_df[redshift_df['PLATE-IFU'].asstr() == plate_ifu]

        if matching_rows.empty:
            print(f"No redshift found for {plate_ifu}. Skipping.")
            continue

        z = matching_rows['REDSHIFT'].iloc[0]
        print(f"Processing file: {filename} (PLATE-IFU: {plate_ifu}, z: {z})")

        # Process this FITS cube
        try:
            # Open FITS file and extract central spectrum
            HDU = fits.open(fits_file)
            flux = HDU[1].data
            wavedata_raw = HDU[6].data * u.AA

            # Get central spaxel coordinates
            x = int(HDU[1].header['CRPIX1'])
            y = int(HDU[1].header['CRPIX2'])
            middleflux = flux[:, x, y]

            # Create output directory
            galaxy_plots_dir = f"{PLOTS_DIR}/{plate_ifu}"
            if not disable_plots:
                os.makedirs(galaxy_plots_dir, exist_ok=True)

            # Measure equivalent widths for all lines
            all_eqw = []
            for i, (line_name, line_wavelength, guess, is_absorption) in enumerate(zip(line_names, lines_to_measure, width_guess, line_absorption)):
                try:
                    plot_full_spectrum = (i == 0) and not disable_plots

                    ew = measure_ew(wavedata_raw, middleflux * u.dimensionless_unscaled,
                                  z, line_wavelength, width_guess=guess,
                                  absorption=is_absorption,
                                  plot=not disable_plots,
                                  plot_full=plot_full_spectrum,
                                  save_dir=galaxy_plots_dir)

                    all_eqw.append(ew.value if hasattr(ew, 'value') else np.nan)
                except Exception as e:
                    print(f"Error measuring {line_name} for {plate_ifu}: {e}")
                    all_eqw.append(np.nan)

            HDU.close()

            # Create results row for this galaxy
            row_dict = {'PLATE-IFU': plate_ifu, 'z': z}
            for name, value in zip(line_names, all_eqw):
                row_dict[name] = value

            # Add to results DataFrame
            results_df = pd.concat([results_df, pd.DataFrame([row_dict])], ignore_index=True)
            measured_count += 1

            # Periodic intermediate saves
            if measured_count % 5 == 0:
                results_df.to_csv(output_csv_path, index=False)
                print(f"Intermediate results saved to {output_csv_path} with {len(results_df)} galaxies")

        except Exception as e:
            print(f"Error processing {plate_ifu}: {e}")
            continue

    # Final save
    if not results_df.empty:
        results_df.to_csv(output_csv_path, index=False)
        print(f"Processing complete. Results saved to {output_csv_path} with {len(results_df)} galaxies")
    else:
        print("No galaxies were successfully processed.")

    return results_df

# COMMAND LINE INTERFACE --------------------------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Measure equivalent widths of spectral lines.')
    parser.add_argument('--preset', choices=list(PRESETS.keys()),
                       help=f'Use a preset configuration: {", ".join(PRESETS.keys())}')
    parser.add_argument('--input-type', choices=['json', 'fits'],
                       help='Type of input files to process')
    parser.add_argument('--folder-path',
                       help='Path to folder containing input files')
    parser.add_argument('--redshift-csv',
                       help='Path to CSV file with redshifts')
    parser.add_argument('--output-csv',
                       help='Path to save output CSV file')
    parser.add_argument('--continuum-method', choices=['improved', 'segmented', 'robust', 'median'],
                       help='Method for continuum fitting')
    parser.add_argument('--plots-dir',
                       help='Directory to save plots')
    parser.add_argument('--no-plots', action='store_true',
                       help='Disable plotting to generate only CSV results')

    args = parser.parse_args()

    # Configure based on preset or individual arguments
    if args.preset:
        preset = PRESETS[args.preset]

        input_type = args.input_type if args.input_type is not None else preset["input_type"]
        folder_path = args.folder_path if args.folder_path is not None else preset["folder_path"]
        redshift_csv = args.redshift_csv if args.redshift_csv is not None else preset["redshift_csv"]
        continuum_method = args.continuum_method if args.continuum_method is not None else preset["continuum_method"]

        # Handle output CSV naming with continuum method
        if args.output_csv is not None:
            output_csv = args.output_csv
        else:
            base_name = os.path.splitext(preset["output_csv"])[0]
            output_csv = f"{base_name}_{continuum_method}.csv"

        # Handle plots directory naming
        if args.plots_dir is not None:
            plots_dir = args.plots_dir
        else:
            plots_dir = preset["plots_dir"]
            if "_" + continuum_method not in plots_dir:
                plots_dir = f"{plots_dir}_{continuum_method}"

    else:
        # Require all essential arguments if no preset specified
        if args.input_type is None or args.folder_path is None or args.redshift_csv is None or args.output_csv is None:
            parser.error("Must specify either --preset or all of: --input-type, --folder-path, --redshift-csv, --output-csv")

        input_type = args.input_type
        folder_path = args.folder_path
        redshift_csv = args.redshift_csv
        continuum_method = args.continuum_method if args.continuum_method is not None else "improved"

        # Handle output CSV naming
        if args.output_csv is not None:
            base_name = os.path.splitext(args.output_csv)[0]
            if "_" + continuum_method not in base_name:
                output_csv = f"{base_name}_{continuum_method}.csv"
            else:
                output_csv = args.output_csv
        else:
            output_csv = f"EW_results_{continuum_method}.csv"

        # Handle plots directory naming
        plots_dir = args.plots_dir if args.plots_dir is not None else "EW_Measurements_Combined"
        if "_" + continuum_method not in plots_dir:
            plots_dir = f"{plots_dir}_{continuum_method}"

    # Update global configuration
    CONTINUUM_METHOD = continuum_method
    PLOTS_DIR = plots_dir

    print(f"\n--- EW Measurement Configuration ---")
    print(f"Input Type: {input_type}")
    print(f"Folder Path: {folder_path}")
    print(f"Redshift CSV: {redshift_csv}")
    print(f"Output CSV: {output_csv}")
    print(f"Continuum Method: {continuum_method}")
    print(f"Plots Directory: {plots_dir}")
    print(f"Plotting: {'Disabled' if args.no_plots else 'Enabled'}")
    print(f"----------------------------------\n")

    # Execute processing based on input file type
    if input_type == 'json':
        process_json_files(folder_path, redshift_csv, output_csv, args.no_plots)
    else:  # fits
        process_fits_files(folder_path, redshift_csv, output_csv, args.no_plots)
