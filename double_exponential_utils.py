#!/usr/bin/env python
"""
Utility functions for working with double exponential quintessence in CAMB and cobaya
"""
import numpy as np
from typing import Dict, List, Tuple

def create_double_exp_params(V0: float = 1e-120, 
                            alpha: float = 0.1, 
                            n: float = 1.0, 
                            V1: float = 5e-121, 
                            beta: float = 0.05,
                            output_file: str = None) -> Dict:
    """
    Create parameters for the double exponential potential:
    V(phi) = V0 * exp(alpha * phi^n) + V1 * exp(beta * phi)
    
    :param V0: amplitude of first exponential term
    :param alpha: coefficient in first exponential
    :param n: power in first exponential
    :param V1: amplitude of second exponential term
    :param beta: coefficient in second exponential
    :param output_file: optional filename to output phi(a) evolution
    :return: dictionary of parameters for set_dark_energy
    """
    params = {
        'dark_energy_model': 'early',  # Uses the EarlyQuintessence class (double exponential only now)
        'V0': V0,
        'alpha': alpha,
        'n': n,
        'V1': V1,
        'beta': beta
    }
    
    if output_file:
        params['output_background_phi'] = True
        params['output_background_phi_filename'] = output_file
    
    return params

def get_cobaya_params() -> Dict:
    """
    Get default parameters for cobaya sampling of double exponential potential
    
    :return: dictionary of parameters for cobaya
    """
    # Create parameters with sensible priors for the double exponential model
    params = {
        # Standard cosmological parameters
        'ombh2': {'prior': {'min': 0.01, 'max': 0.03}, 
                  'ref': {'dist': 'norm', 'loc': 0.0224, 'scale': 0.0002}, 
                  'latex': '\\Omega_b h^2'},
                  
        'omch2': {'prior': {'min': 0.05, 'max': 0.25}, 
                 'ref': {'dist': 'norm', 'loc': 0.12, 'scale': 0.002}, 
                 'latex': '\\Omega_c h^2'},
                 
        'H0': {'prior': {'min': 60, 'max': 80}, 
              'ref': {'dist': 'norm', 'loc': 67.5, 'scale': 1.0}, 
              'latex': 'H_0'},
              
        # Double exponential model parameters
        'V0': {'prior': {'min': 0.1e-120, 'max': 10e-120}, 
              'ref': {'dist': 'norm', 'loc': 1e-120, 'scale': 1e-121}, 
              'latex': 'V_0'},
              
        'alpha': {'prior': {'min': 0.01, 'max': 0.5}, 
                 'ref': {'dist': 'norm', 'loc': 0.1, 'scale': 0.05}, 
                 'latex': '\\alpha'},
                 
        'n': {'prior': {'min': 0.5, 'max': 2.0}, 
             'ref': {'dist': 'norm', 'loc': 1.0, 'scale': 0.2}, 
             'latex': 'n'},
             
        'V1': {'prior': {'min': 0.1e-121, 'max': 10e-121}, 
              'ref': {'dist': 'norm', 'loc': 5e-121, 'scale': 1e-121}, 
              'latex': 'V_1'},
              
        'beta': {'prior': {'min': 0.01, 'max': 0.2}, 
                'ref': {'dist': 'norm', 'loc': 0.05, 'scale': 0.02}, 
                'latex': '\\beta'}
    }
    
    return params

def create_cobaya_yaml(output_file: str = 'double_exp_cobaya.yaml') -> None:
    """
    Create a cobaya YAML file for sampling the double exponential model
    
    :param output_file: output filename for the YAML configuration
    """
    import yaml
    
    # Get parameter priors
    params = get_cobaya_params()
    
    # Create cobaya configuration
    config = {
        'params': params,
        'likelihood': {
            # Example likelihoods - replace with ones you want to use
            'planck_2018_lowl.TT': None,
            'planck_2018_lowl.EE': None,
            'planck_2018_highl_plik.TTTEEE': None
        },
        'theory': {
            'camb': {
                'extra_args': {
                    'dark_energy_model': 'early',
                    'lens_potential_accuracy': 1
                }
            }
        },
        'sampler': {
            'mcmc': {
                'Rminus1_stop': 0.01,
                'max_samples': 50000
            }
        },
        'output': 'chains/double_exp'
    }
    
    # Write to file
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    print(f"Cobaya configuration saved to {output_file}")

if __name__ == "__main__":
    # Example usage
    create_cobaya_yaml()
    print("\nExample parameters for CAMB:")
    params = create_double_exp_params(V0=1e-120, alpha=0.1, n=1.0, V1=5e-121, beta=0.05)
    for k, v in params.items():
        print(f"  {k}: {v}") 