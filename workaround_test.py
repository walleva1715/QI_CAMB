#!/usr/bin/env python3
"""
Test the early Quintessence dark energy model via the CAMB python wrapper
with a workaround for the search_for_initialphi issue.
"""
import camb
from camb.model import CAMBparams
import numpy as np
import os

def main():
    # Initialize CAMB parameters from existing ini file known to work
    ini_path = os.path.expanduser("~/QI_CAMB/fortran/Modified_files/params_quintessence.ini")
    params = camb.read_ini(ini_path)
    
    # Get a reference to the dark energy model
    de = params.DarkEnergy
    print(f"Loaded dark energy model: {type(de).__name__}")
    print(f"Potential type: {de.potential_type}")
    print(f"Potential params: {[de.potentialparams[i] for i in range(5)]}")
    
    # Calculate background cosmology
    # Note: we'll need to verify params since they were read from the .ini file
    print(f"Cosmological parameters: H0={params.H0}, ombh2={params.ombh2}, omch2={params.omch2}")
    
    # Calculate results
    results = camb.get_background(params)
    
    # Print some results
    #print(f"Age of the universe: {results.age_of_universe:.2f} Gyr")
    print(f"Dark energy density today: {results.omega_de:.5f}")
    
    # Get the scale factor array
    a_array = np.logspace(-4, 0, 100)
    # Get the corresponding redshifts
    z_array = 1/a_array - 1
    
    # Calculate Hubble rate at each redshift
    h_of_z = results.h_of_z(z_array)
    
    print(f"H(z=0): {h_of_z[-1]:.5f}")
    print(f"H(z=1): {h_of_z[np.argmin(np.abs(z_array - 1.0))]:.5f}")
    
    print("Test completed successfully!")

if __name__ == "__main__":
    main() 