"""
Example of using Quintessential Inflation with Cobaya for parameter inference.

This script demonstrates how to set up a custom Cobaya Theory provider
that interfaces with CAMB's Quintessential Inflation model.
"""

import numpy as np
import camb
from camb import model, dark_energy

# Define a custom Theory provider class for Cobaya
class QuintessentialTheory:
    """
    Custom theory provider for Cobaya to interface with CAMB's
    Quintessential Inflation model.
    
    For actual use, you would import Cobaya's theory classes.
    """
    
    def initialize(self):
        """Set up the class when Cobaya initializes it"""
        # We'll keep track of current parameters and results for efficiency
        self.current_params = None
        self.current_results = None
    
    def get_requirements(self):
        """Tell Cobaya what parameters we need"""
        # Basic cosmology
        return {"H0", "ombh2", "omch2", 
                # Quintessence parameters
                "QI_potential_type", 
                "QI_param1", "QI_param2", "QI_param3"}
    
    def calculate(self, state, want_derived=True, **params_values):
        """Calculate results for a given set of parameters"""
        # Only recalculate if parameters have changed
        H0 = params_values["H0"]
        ombh2 = params_values["ombh2"]
        omch2 = params_values["omch2"]
        pot_type = params_values["QI_potential_type"]
        pot_params = [
            params_values["QI_param1"],
            params_values["QI_param2"],
            params_values["QI_param3"],
            0.0,  # Fixed parameter 4
            0.0,  # Fixed parameter 5
        ]
        
        # New parameter set - needs calculation
        cambparams = camb.CAMBparams()
        cambparams.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
        
        # Set up quintessence parameters directly with set_dark_energy
        cambparams.set_dark_energy(
            dark_energy_model='early',
            n=3.0,  # Fixed exponent
            potential_type=pot_type,
            potentialparams=pot_params
        )
        
        # Run CAMB
        self.current_results = camb.get_results(cambparams)
        # Store parameters we used
        self.current_params = params_values.copy()
    
    def get_Omega(self, component, **params_values):
        """Get Omega values for different components"""
        if self.current_results is None:
            raise ValueError("Must calculate() first")
        return self.current_results.get_Omega(component)
    
    def get_derived_parameters(self, derived=None, **params_values):
        """Return derived parameters that Cobaya will report"""
        if self.current_results is None:
            raise ValueError("Must calculate() first")
        
        # Return derived parameters
        derived_params = {
            "Omega_de": self.current_results.get_Omega("de"),
            "Omega_m": self.current_results.get_Omega("m"),
            "sigma8": self.current_results.get_sigma8()
        }
        return derived_params


# Example Cobaya-like interface
if __name__ == "__main__":
    # This would be replaced with real Cobaya code
    print("Example of Quintessential Inflation Theory for Cobaya")
    
    # Create a Theory instance
    theory = QuintessentialTheory()
    theory.initialize()
    
    # Create example parameter set
    test_params = {
        "H0": 67.5,
        "ombh2": 0.022,
        "omch2": 0.122,
        "QI_potential_type": 7,
        "QI_param1": 7e-121,
        "QI_param2": 2.0,
        "QI_param3": 2.0
    }
    
    # Run calculation
    print("\nRunning calculation with parameters:")
    for k, v in test_params.items():
        print(f"  {k}: {v}")
    
    # Create dummy state object for Cobaya compatibility
    class DummyState:
        pass
    
    # Call calculate 
    theory.calculate(DummyState(), **test_params)
    
    # Get derived parameters
    derived = theory.get_derived_parameters(**test_params)
    
    # Print results
    print("\nDerived parameters:")
    for k, v in derived.items():
        print(f"  {k}: {v:.5f}")
    
    print("\nExample completed successfully!") 