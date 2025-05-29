import camb
import numpy as np
from camb import model, initialpower, dark_energy

print("CAMB Quintessential Inflation Test")

# EXAMPLE 1: Direct use of set_dark_energy with parameters
print("\nEXAMPLE 1: Using set_dark_energy with parameters")
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, omk=0)

# Set up the quintessence model with parameters
pars.set_dark_energy(
    dark_energy_model='early',        # Use the EarlyQuintessence model
    potential_type=7,                 # General exponential potential 
    potentialparams=[1.0e-5, 2.0, 2.0, 0.0, 0.0],  # Use more reasonable value for V_0
    output_background_phi=True,
    output_background_phi_filename="phi_evolution.dat"
)

# Print parameters to confirm they're set correctly
print("Potential type:", pars.DarkEnergy.potential_type)
print("Potential params:", [pars.DarkEnergy.potentialparams[i] for i in range(5)])
print("search_for_initialphi:", pars.DarkEnergy.search_for_initialphi)

# Run CAMB with these parameters
try:
    results = camb.get_results(pars)
    print("✓ CAMB calculation successful!")
    print(f"Omega_de = {results.get_Omega('de'):.5f}")
    print(f"H(z=0) = {results.hubble_parameter(0):.8f}")
    print(f"H(z=1) = {results.hubble_parameter(1):.8f}")
except Exception as e:
    print(f"✗ CAMB calculation failed: {e}")


# EXAMPLE 2: Creating the model first, then passing to set_dark_energy
print("\nEXAMPLE 2: Creating model instance first")
pars2 = camb.CAMBparams()
pars2.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, omk=0)

# Create the quintessence model directly
qi_model = dark_energy.EarlyQuintessence()

# Set potential type and params (explicitly set each parameter)
qi_model.potential_type = 7
qi_model.potentialparams[0] = 1.0e-5  # More reasonable V_0 value
qi_model.potentialparams[1] = 2.0
qi_model.potentialparams[2] = 2.0
qi_model.potentialparams[3] = 0.0
qi_model.potentialparams[4] = 0.0

# Set other parameters
qi_model.output_background_phi = True
qi_model.output_background_phi_filename = b"phi_evolution2.dat"  # Note: needs to be bytes
qi_model.search_for_initialphi = False

# Verify parameter values
print("Potential type:", qi_model.potential_type)
print("Potential params:", [qi_model.potentialparams[i] for i in range(5)])
print("search_for_initialphi:", qi_model.search_for_initialphi)

# Set the dark energy model
pars2.set_dark_energy(dark_energy_model=qi_model)

# Run CAMB with these parameters
try:
    results2 = camb.get_results(pars2)
    print("✓ CAMB calculation successful!")
    print(f"Omega_de = {results2.get_Omega('de'):.5f}")
    print(f"H(z=0) = {results2.hubble_parameter(0):.8f}")
    print(f"H(z=1) = {results2.hubble_parameter(1):.8f}")
except Exception as e:
    print(f"✗ CAMB calculation failed: {e}")

print("\nTest completed!") 