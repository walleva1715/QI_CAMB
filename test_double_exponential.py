#!/usr/bin/env python
"""
Test script for double exponential quintessence model
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Make sure we can import from the current directory
sys.path.insert(0, os.path.abspath('.'))

import camb
from camb import model

def test_double_exponential():
    """Test the double exponential quintessence potential"""
    
    # Create CAMB parameters
    pars = camb.CAMBparams()
    
    # Set basic cosmology
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12)
    
    # Configure double exponential quintessence
    # These parameters are chosen to be conservative and enhance stability
    pars.set_dark_energy(
        dark_energy_model='early',  # Uses EarlyQuintessence class (always double exponential now)
        potentialparams=[9e-120, 0.1, 1.0, 5e-121, 0.05],  # [V0, alpha, n, V1, beta]
        output_background_phi=True,
        output_background_phi_filename="quintessence_background.dat"
    )
    
    # Set up for CMB calculation
    pars.set_for_lmax(2500)
    
    # Run CAMB
    print("Running CAMB with double exponential quintessence...")
    results = camb.get_results(pars)
    
    # Check dark energy properties
    background = results.get_background()
    w_z0 = background.get_dark_energy_w_de(0)
    print(f"Dark energy equation of state at z=0: w(z=0) = {w_z0:.6f}")
    
    # Get Omega_de
    omega_de = results.Params.Omega_de
    print(f"Dark energy density parameter: Omega_de = {omega_de:.6f}")
    
    # Plot the equation of state evolution
    zs = np.logspace(-3, 2, 100)
    w_z = np.array([background.get_dark_energy_w_de(z) for z in zs])
    
    plt.figure(figsize=(10, 6))
    plt.semilogx(zs, w_z)
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.xlabel("Redshift z", fontsize=12)
    plt.ylabel("w(z)", fontsize=12)
    plt.title("Dark Energy Equation of State - Double Exponential Potential", fontsize=14)
    plt.savefig("w_evolution.png")
    print("Equation of state evolution saved to w_evolution.png")
    
    # Check to see if the background file was created
    if os.path.exists("quintessence_background.dat"):
        print("Background evolution file was successfully created")
        # Plot the background evolution (first few columns are a, phi, phidot)
        try:
            bg_data = np.loadtxt("quintessence_background.dat", skiprows=1)
            plt.figure(figsize=(12, 8))
            
            # Top panel: phi(a)
            plt.subplot(2, 1, 1)
            plt.plot(bg_data[:, 0], bg_data[:, 1])
            plt.xscale('log')
            plt.grid(True)
            plt.ylabel(r'$\phi$', fontsize=12)
            plt.title('Scalar Field Evolution', fontsize=14)
            
            # Bottom panel: w(a)
            plt.subplot(2, 1, 2)
            plt.plot(bg_data[:, 0], bg_data[:, 4])
            plt.xscale('log')
            plt.grid(True)
            plt.xlabel('Scale Factor a', fontsize=12)
            plt.ylabel('w(a)', fontsize=12)
            
            plt.tight_layout()
            plt.savefig("phi_evolution.png")
            print("Scalar field evolution saved to phi_evolution.png")
        except Exception as e:
            print(f"Could not plot background file: {e}")
    
    # Return results for further analysis if needed
    return results

if __name__ == "__main__":
    try:
        results = test_double_exponential()
        print("\nTest completed successfully!")
    except Exception as e:
        print(f"Test failed with error: {str(e)}") 