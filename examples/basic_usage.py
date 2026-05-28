"""
Basic MEWS Usage Example
========================
Demonstrates simple equivalent width measurement for a single galaxy spectrum.

Author: Olivia A. Greene, PhD
Contact: oliviaallegragreene@gmail.com
"""

import sys
sys.path.append('..')  # Add parent directory to path

import numpy as np
from astropy import units as u

# Import MEWS functions
# NOTE: Update this import based on your actual module structure
# from MEWS_Measuring_Equivalent_Width_in_Spectra import measure_ew

def basic_example():
    """
    Demonstrate basic EW measurement workflow
    """
    
    print("="*60)
    print("BASIC MEWS USAGE EXAMPLE")
    print("="*60)
    
    print("\nThis example shows how to:")
    print("  1. Load spectrum data")
    print("  2. Measure equivalent width for a spectral line")
    print("  3. Generate diagnostic plots")
    
    # ========================================================================
    # STEP 1: Prepare Your Data
    # ========================================================================
    
    print("\n" + "="*60)
    print("STEP 1: Load Spectrum Data")
    print("="*60)
    
    # Example wavelength array (replace with your actual data)
    wavelength = np.linspace(3500, 10000, 5000) * u.AA
    
    # Example flux array (replace with your actual data)
    # For real data, load from FITS or JSON file
    flux = np.random.randn(5000) + 100  
    flux_units = u.dimensionless_unscaled
    
    # Your galaxy's redshift
    redshift = 0.05
    
    print(f"  Wavelength range: {wavelength[0]:.1f} - {wavelength[-1]:.1f}")
    print(f"  Number of pixels: {len(wavelength)}")
    print(f"  Galaxy redshift: {redshift:.4f}")
    
    # ========================================================================
    # STEP 2: Define Spectral Line to Measure
    # ========================================================================
    
    print("\n" + "="*60)
    print("STEP 2: Select Spectral Line")
    print("="*60)
    
    # Balmer lines (common for E+A galaxies)
    H_alpha = 6564.61 * u.AA    # Absorption in E+A galaxies
    H_beta = 4862.68 * u.AA     # Absorption
    H_gamma = 4341.68 * u.AA    # Absorption
    H_delta = 4102.89 * u.AA    # Absorption (strong E+A diagnostic)
    
    # For this example, measure H-delta
    line_wavelength = H_delta
    width_guess = 3.0  # Initial width estimate in Angstroms
    
    print(f"  Measuring: H-delta")
    print(f"  Rest wavelength: {line_wavelength:.2f}")
    print(f"  Observed wavelength (z={redshift}): {line_wavelength * (1 + redshift):.2f}")
    
    # ========================================================================
    # STEP 3: Measure Equivalent Width
    # ========================================================================
    
    print("\n" + "="*60)
    print("STEP 3: Measure Equivalent Width")
    print("="*60)
    
    # Uncomment and modify once you import the actual function:
    """
    ew_result = measure_ew(
        wavelength=wavelength,
        flux=flux * flux_units,
        z=redshift,
        line_wavelength=line_wavelength,
        width_guess=width_guess,
        absorption=True,  # True for absorption lines (H-delta, H-gamma, etc.)
        plot=True,  # Generate diagnostic plots
        save_dir="output_plots",
        continuum_method='improved'  # Use validated method
    )
    
    print(f"\n  Measured EW: {ew_result:.2f} Angstroms")
    print(f"  Diagnostic plot saved to: output_plots/")
    """
    
    print("\n  [Replace example data with your actual spectrum to run measurement]")
    
    # ========================================================================
    # STEP 4: Batch Processing (Multiple Galaxies)
    # ========================================================================
    
    print("\n" + "="*60)
    print("STEP 4: Batch Processing (Optional)")
    print("="*60)
    
    print("\n  For multiple galaxies, use the CLI interface:")
    print("\n  python MEWS_Measuring_Equivalent_Width_in_Spectra.py \\")
    print("    --preset example_json \\")
    print("    --no-plots")
    
    print("\n" + "="*60)
    print("EXAMPLE COMPLETE")
    print("="*60)
    print("\nNext steps:")
    print("  1. Replace example data with your actual spectrum")
    print("  2. Uncomment the measure_ew() function call")
    print("  3. Run this script to see EW measurements")
    print("  4. Check output_plots/ for diagnostic figures")
    print("\nQuestions? Contact: oliviaallegragreene@gmail.com")
    print("="*60)


if __name__ == "__main__":
    basic_example()
