"""
Quintessence Inflation Theory implementation for Cobaya with custom CAMB interface.

This file implements a Theory provider for Cobaya that computes quintessence inflation
with the potential V(φ) = V0·exp[-γ φ^n_exp], supporting both parameter inference and
prediction of inflationary observables (n_s, r).
"""

from cobaya.theory import Theory
import numpy as np
import camb
from camb import model, dark_energy
import os

class QuintessenceInflationTheory(Theory):
    """
    Quintessence Inflation Theory with potential V(φ) = V0·exp[-γ φ^n_exp].
    
    This theory provider computes:
    1. Analytical slow-roll predictions for inflationary observables
    2. CMB power spectra via customized CAMB integration
    
    Analytical slow-roll derivation:
    
    1) End of inflation: ε(φ_end)=1 ⇒ φ_end = [√2/(γ·n_exp)]^(1/(n_exp−1)).
    
    2) Number of e-folds from φ_end to φ_star:
       For n_exp≠2: N = (1/(γ·n_exp·(n_exp−2)))·(φ_star^(−(n_exp−2)) − φ_end^(−(n_exp−2)))
         ⇒ φ_star = [φ_end^(−(n_exp−2)) − γ·n_exp·(n_exp−2)·N]^(−1/(n_exp−2)).
       For n_exp=2: N = γ·(φ_star^2 − φ_end^2)/4 ⇒ φ_star = sqrt(φ_end^2 + 4N/γ).
    
    3) Slow-roll parameters at φ_star:
       ε = ½ [γ·n_exp·φ_star^(n_exp−1)]^2,
       η = γ·n_exp·(n_exp−1)·φ_star^(n_exp−2) − [γ·n_exp·φ_star^(n_exp−1)]^2.
    
    4) Observables:
       n_s = 1 − 6ε + 2η = 1 − (1/N)[2n_exp−2 + γ·n_exp·(1/(γ·n_exp·(n_exp−2)·N))^(n_exp/(n_exp−2))],
       r   = 16ε   = 8·n_exp^2·γ^2·(1/(γ·n_exp·(n_exp−2)·N))^((2n_exp−2)/(n_exp−2)).
    """
    
    # Define model parameters
    params = {
        # Free model parameters
        "N_efolds": None,
        "gamma": None,
        "n_exp": None,
        
        # Common cosmological parameters that will be passed to CAMB
        "H0": None,         # Hubble constant
        "ombh2": None,      # Baryon density
        "omch2": None,      # Cold dark matter density
        "tau": None,        # Optical depth
        "As": None,         # Scalar amplitude
        
        # Derived inflationary observables
        "n_s": {"derived": True},
        "r": {"derived": True},
        
        # Optional additional derived parameters
        "Omega_de": {"derived": True},
        "Omega_m": {"derived": True},
        "sigma8": {"derived": True}
    }

    def initialize(self):
        """Initialize the theory code when Cobaya starts"""
        self.current_params = None
        self.current_results = None

    def get_requirements(self):
        """Tell Cobaya what we need from the sampler"""
        return {}

    def get_can_provide(self):
        """Tell Cobaya what this theory can compute"""
        return ["Cl", "H0", "n_s", "r", "Omega_de", "Omega_m", "sigma8"]

    def calculate(self, state, want_derived=True, **params_values):
        """
        Main calculation method called by Cobaya.
        
        Computes n_s and r analytically from the potential parameters,
        then runs CAMB to get power spectra and other derived parameters.
        """
        # 1. Extract parameters
        N = params_values["N_efolds"]
        gamma = params_values["gamma"]
        n_exp = params_values["n_exp"]
        
        # 2. Compute n_s and r from analytical slow-roll
        n_s, r = self.compute_ns_r_from_model(N, gamma, n_exp)
        # Check if compute_ns_r_from_model indicated an error or returned unphysical values
        unphysical_ns = not (0.8 < n_s < 1.2) # Typical range for n_s
        unphysical_r = not (0 <= r < 1.0)  # Typical range for r (r must be non-negative)
        if np.isinf(n_s) or np.isinf(r) or np.isnan(n_s) or np.isnan(r) or unphysical_ns or unphysical_r:
            self.log.debug(
                "Unphysical n_s/r from model: N_efolds=%s, gamma=%s, n_exp=%s -> n_s=%s, r=%s. Rejecting point.",
                N, gamma, n_exp, n_s, r)
            state["logp"] = -np.inf
            return False # Indicate failure to Cobaya
        state["n_s"] = n_s
        state["r"] = r
        
        # 3. Set up CAMB parameters
        cambparams = camb.CAMBparams()
        
        # Base cosmology
        H0 = params_values.get("H0", 67.5)
        ombh2 = params_values.get("ombh2", 0.022)
        omch2 = params_values.get("omch2", 0.120)
        tau = params_values.get("tau", 0.06)
        As = params_values.get("As", 2.1e-9)
        
        cambparams.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, tau=tau)
        
        # 4. Set up initial power spectrum with n_s and r from our model
        cambparams.InitPower.set_params(As=As, ns=n_s, r=r)
        
        # Ensure tensor modes are computed when r>0
        if r > 0.0:
            cambparams.WantTensors = True
        
        # 5. Configure the quintessence potential
        # For the exponential potential V(φ) = V0·exp[-γ φ^n_exp], use potential_type 7
        V0 = 9e-121  # Placeholder V0 value - ensures quintessence energy density is low during inflation
        potentialparams = [V0, gamma, n_exp, 0.0, 0.0]
        
        # Set dark energy model
        cambparams.set_dark_energy(
            dark_energy_model='early',
            potential_type=7,  # Generalized exponential
            potentialparams=potentialparams
        )
        
        # 6. Set accuracy (disable non-linear Halofit to prevent HMCode integration timeout)
        cambparams.set_for_lmax(lmax=2500, lens_potential_accuracy=1, nonlinear=False)
        
        # 7. Request matter transfer for sigma8 (ensure WantTransfer=True)
        cambparams.set_matter_power(redshifts=[0.0], kmax=2.0, nonlinear=False)
        # 8. Run CAMB
        self.current_results = camb.get_results(cambparams)
        
        # 9. Store derived parameters
        if want_derived:
            state["derived"] = {
                "n_s": n_s,
                "r": r,
                "Omega_de": self.current_results.get_Omega("de"),
                "Omega_m": self.current_results.get_Omega("cdm"),
                "sigma8": self.current_results.get_sigma8()
            }
        
        # 10. Store power spectra
        powers = self.current_results.get_cmb_power_spectra()
        state["Cl"] = powers["total"]
        
        # Store current parameters for caching
        self.current_params = params_values.copy()

    def get_Cl(self, **kwargs):
        """Return CMB power spectra as dict mapping spectral names to arrays for Cobaya."""
        # Retrieve all power spectra from CAMB results
        powers = self.current_results.get_cmb_power_spectra()
        total = powers["total"]
        cl_dict = {}
        # Map CAMB output columns: TT, EE, BB, TE
        cl_dict["tt"] = total[:, 0]
        cl_dict["ee"] = total[:, 1]
        cl_dict["bb"] = total[:, 2]
        cl_dict["te"] = total[:, 3]
        # Include lensing potential spectrum if available
        if "lens_potential" in powers:
            lens = powers["lens_potential"]
            # lens may be 1D or 2D array
            cl_dict["pp"] = lens if lens.ndim == 1 else lens[:, 0]
        # Default zero arrays for tb and eb
        cl_dict["tb"] = np.zeros_like(cl_dict["tt"])
        cl_dict["eb"] = np.zeros_like(cl_dict["ee"])
        return cl_dict

    def compute_ns_r_from_model(self, N, gamma, n_exp):
        """
        Slow-roll prediction for V(φ)=V0*exp[-γ φ^n_exp] in M_pl=1 units.
        Returns (n_s, r).
        """
    
        # 1) φ_end from ε=1: 0.5*(γ n φ_end^(n-1))^2 = 1
        phi_end = (np.sqrt(2)/(abs(gamma)*n_exp))**(1.0/(n_exp-1))
    
        # 2) invert e-fold integral
        if abs(n_exp - 2.0) > 1e-10:
            # General n ≠ 2
            # φ_*^(2-n) = φ_end^(2-n) - γ n (2-n) N
            power = 2.0 - n_exp
            inv = phi_end**power - gamma * n_exp * power * N
            phi_star = inv**(1.0/power)
        else:
            # n = 2  ⇒  φ_* = φ_end * exp(-2 γ N)
            phi_star = phi_end * np.exp(-2.0 * gamma * N)
    
        # 3) slow-roll parameters
        # ε = 0.5*(γ n φ^(n-1))^2
        epsilon = 0.5 * (gamma * n_exp * phi_star**(n_exp-1))**2
    
        # η = g'' + (g')^2  with g' = -γn φ^(n-1), g'' = -γn(n-1)φ^(n-2)
        term_g2 = (gamma * n_exp * phi_star**(n_exp-1))**2
        term_gpp = -gamma * n_exp * (n_exp-1) * phi_star**(n_exp-2)
        eta = term_gpp + term_g2
    
        # 4) observables
        n_s = 1.0 - 6.0*epsilon + 2.0*eta
        r   = 16.0*epsilon
    
        return n_s, r


# Example usage
if __name__ == "__main__":
    print("Example of Quintessence Inflation Theory for Cobaya")
    
    # Create a Theory instance
    theory = QuintessenceInflationTheory()
    theory.initialize()
    
    # Create example parameter set
    test_params = {
        "H0": 67.5,
        "ombh2": 0.022,
        "omch2": 0.122,
        "tau": 0.06,
        "As": 2.1e-9,
        "N_efolds": 50.0,
        "gamma": 0.1,
        "n_exp": 2.0
    }
    
    # Run calculation
    print("\nRunning calculation with parameters:")
    for k, v in test_params.items():
        print(f"  {k}: {v}")
    
    # Create dummy state object for Cobaya compatibility
    class DummyState(dict):
        pass
    
    # Call calculate 
    state = DummyState()
    theory.calculate(state, **test_params)
    
    # Print results
    print("\nInflationary observables:")
    print(f"  n_s: {state['n_s']:.6f}")
    print(f"  r: {state['r']:.6f}")
    
    if "derived" in state:
        print("\nDerived parameters:")
        for k, v in state["derived"].items():
            if isinstance(v, float):
                print(f"  {k}: {v:.6f}")
    
    print("\nExample completed successfully!") 