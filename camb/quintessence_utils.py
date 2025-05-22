"""
Helper utilities for working with quintessence models in CAMB and cobaya
"""

import numpy as np
from typing import List, Dict, Union, Optional

def create_quintessence_params(potential_type: int, parameters: List[float], 
                              output_file: Optional[str] = None) -> Dict:
    """
    Helper to create quintessence parameters for use with CAMB and cobaya
    
    :param potential_type: Integer (0-8) selecting the potential:
        0: Early quintessence axion-like (1-cos)^n potential
        1: Harmonic potential (thawing) m²φ²/2
        2: Inverse potential (freezing) M^5/φ
        3: Cubic potential m*φ^3/3
        4: Inverse square potential M/φ²
        5: General power law potential m*φ^n
        6: Hyperbolic cosine well V₀*cosh(β*(φ/Mpl)^u)
        7: Generalized exponential V₀*exp(α*φ^n)
        8: Double exponential V₀*exp(α*φ^n) + V₁*exp(β*φ)
        
    :param parameters: List of parameters for the chosen potential:
        - For type 1: [m]
        - For type 2: [m]
        - For type 3: [m]
        - For type 4: [m]
        - For type 5: [m, n]
        - For type 6: [V₀, β, u]
        - For type 7: [V₀, α, n]
        - For type 8: [V₀, α, n, V₁, β]
        
    :param output_file: Optional filename to output background evolution
    :return: Dictionary for use with CAMBparams.set_dark_energy
    """
    params = {
        'dark_energy_model': 'early',
        'potential_type': potential_type,
        'potentialparams': parameters
    }
    
    if output_file:
        params['output_background_phi'] = True
        params['output_background_phi_filename'] = output_file
        
    return params

def create_cobaya_priors(potential_type: int) -> Dict:
    """
    Create suggested priors for cobaya sampling based on the quintessence potential type.
    These are just suggestions and should be adjusted for your specific analysis.
    
    :param potential_type: Integer (0-8) selecting the potential
    :return: Dictionary of parameter priors for cobaya
    """
    priors = {}
    
    if potential_type == 0:  # Early quintessence axion-like
        priors = {
            'n': {'prior': [1.0, 5.0], 'ref': {'dist': 'norm', 'loc': 3.0, 'scale': 0.5}},
            'f': {'prior': [0.01, 0.1], 'ref': {'dist': 'norm', 'loc': 0.05, 'scale': 0.01}},
            'm': {'prior': [1e-54, 1e-53], 'ref': {'dist': 'norm', 'loc': 5e-54, 'scale': 1e-54}},
            'theta_i': {'prior': [0.1, 3.14], 'ref': {'dist': 'norm', 'loc': 3.0, 'scale': 0.5}}
        }
    elif potential_type == 1:  # Harmonic potential
        priors = {
            'pot_m': {'prior': [1e-120, 1e-110], 'ref': {'dist': 'norm', 'loc': 1e-115, 'scale': 1e-116}}
        }
    elif potential_type == 2:  # Inverse potential
        priors = {
            'pot_m': {'prior': [1e-120, 1e-110], 'ref': {'dist': 'norm', 'loc': 1e-115, 'scale': 1e-116}}
        }
    elif potential_type == 5:  # General power law
        priors = {
            'pot_m': {'prior': [1e-120, 1e-110], 'ref': {'dist': 'norm', 'loc': 1e-115, 'scale': 1e-116}},
            'pot_n': {'prior': [1.0, 4.0], 'ref': {'dist': 'norm', 'loc': 2.0, 'scale': 0.5}}
        }
    elif potential_type == 7:  # Generalized exponential
        priors = {
            'V0': {'prior': [0.5e-120, 5e-120], 'ref': {'dist': 'norm', 'loc': 1e-120, 'scale': 0.5e-120}},
            'alpha': {'prior': [0.5, 2.0], 'ref': {'dist': 'norm', 'loc': 1.0, 'scale': 0.2}},
            'exp_n': {'prior': [1.0, 3.0], 'ref': {'dist': 'norm', 'loc': 2.0, 'scale': 0.5}}
        }
    elif potential_type == 8:  # Double exponential
        priors = {
            'V0': {'prior': [0.5e-120, 5e-120], 'ref': {'dist': 'norm', 'loc': 1e-120, 'scale': 0.5e-120}},
            'alpha': {'prior': [0.5, 2.0], 'ref': {'dist': 'norm', 'loc': 1.0, 'scale': 0.2}},
            'exp_n': {'prior': [1.0, 3.0], 'ref': {'dist': 'norm', 'loc': 2.0, 'scale': 0.5}},
            'V1': {'prior': [0.5e-120, 5e-120], 'ref': {'dist': 'norm', 'loc': 1e-120, 'scale': 0.5e-120}},
            'beta': {'prior': [0.5, 2.0], 'ref': {'dist': 'norm', 'loc': 1.0, 'scale': 0.2}}
        }
    
    return priors

def create_cobaya_example(potential_type: int = 7) -> Dict:
    """
    Create an example cobaya input for the selected quintessence potential
    
    :param potential_type: Integer (0-8) selecting the potential
    :return: Example cobaya configuration dictionary
    """
    priors = create_cobaya_priors(potential_type)
    
    # Basic cobaya configuration
    config = {
        'params': {
            # Standard cosmology parameters
            'ombh2': 0.0224, 
            'omch2': 0.120,
            'H0': 67.5,
            'tau': 0.06,
            'As': 2.1e-9,
            'ns': 0.965,
            # Add quintessence parameters from priors
            **priors
        },
        'likelihood': {
            # Example likelihoods - replace with actual ones
            'planck_2018_lowl.TT': None,
            'planck_2018_lowl.EE': None,
            'planck_2018_highl_plik.TTTEEE': None,
            'planck_2018_lensing.clik': None
        },
        'theory': {
            'camb': {
                # Set base accuracy and configure quintessence
                'extra_args': {
                    'lens_potential_accuracy': 1,
                    'DoLateRadTruncation': False,
                    'dark_energy_model': 'early',
                    'potential_type': potential_type
                }
            }
        },
        'sampler': {
            'mcmc': {
                'Rminus1_stop': 0.01,
                'max_samples': 10000
            }
        }
    }
    
    return config

# Example of how to use these utilities in a cobaya likelihood
def example_cobaya_usage():
    """
    Example showing how to use these utilities with cobaya
    """
    code = """
    # Example cobaya likelihood using quintessence
    from cobaya.likelihood import Likelihood
    from cobaya.theory import Theory
    import camb
    from camb import model, quintessence_utils
    import numpy as np
    
    class QuintessenceLikelihood(Likelihood):
        # Define parameters needed from sampler
        params = {
            "ombh2": None,
            "omch2": None, 
            "H0": None,
            "V0": None,
            "alpha": None,
            "exp_n": None
        }
        
        # Initialize
        def initialize(self):
            self.potential_type = 7  # Generalized exponential
            
        # Calculate likelihood
        def logp(self, **params_values):
            # Get theory predictions
            theory_results = self.provider.get_theory_results()
            
            # Calculate likelihood based on theory_results...
            # (This is just an example)
            
            return 0.0  # Return log-likelihood
            
    class QuintessenceTheory(Theory):
        # Define parameters needed from sampler
        params = {
            "ombh2": None,
            "omch2": None, 
            "H0": None,
            "V0": None,
            "alpha": None,
            "exp_n": None
        }
        
        # Initialize
        def initialize(self):
            self.potential_type = 7  # Generalized exponential
            
        # Get results
        def get_theory_results(self, **params_values):
            # Set up CAMB
            camb_params = camb.CAMBparams()
            
            # Set basic cosmology
            camb_params.set_cosmology(
                H0=params_values['H0'], 
                ombh2=params_values['ombh2'], 
                omch2=params_values['omch2']
            )
            
            # Set quintessence potential
            potentialparams = [
                params_values['V0'],
                params_values['alpha'],
                params_values['exp_n']
            ]
            
            # Use the utility function to configure quintessence
            quint_params = quintessence_utils.create_quintessence_params(
                potential_type=self.potential_type, 
                parameters=potentialparams
            )
            
            # Apply to CAMB parameters
            camb_params.set_dark_energy(**quint_params)
            
            # Set other CAMB parameters
            camb_params.set_for_lmax(2500, lens_potential_accuracy=1)
            
            # Calculate results
            results = camb.get_results(camb_params)
            powers = results.get_cmb_power_spectra()
            
            # Return dictionary of results for likelihood function
            return {
                'powers': powers,
                'background': results.get_background()
            }
    """
    return code 