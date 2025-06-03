from cobaya.theory import Theory
import numpy as np
import camb
from scipy import constants

class QICAMBTheory(Theory):
    """
    Link Cobaya to QI CAMB with analytic inflation+quintessence DE.

    Inputs:
      - N_efolds: float
      - gamma: float
      - n_exp: float
      - H0: float
      - ombh2: float
      - omch2: float
      - tau: float
      - As: float
      - mnu: float
      - YHe: float
      - potentialparam1..5: double-exponential quintessence DE

    Provides:
      Cl, n_s, r, Omega_de, omegam, sigma8, age, rdrag, zrei, Y_p, DH
    """
    params = {
        "N_efolds": None,
        "gamma": None,
        "n_exp": None,
        "H0": None,
        "ombh2": None,
        "omch2": None,
        "tau": None,
        "As": None,
        "mnu": None,
        "YHe": None,
        "potentialparam1": 9.0e-151,
        "potentialparam2": 2,
        "potentialparam3": 4,
        "potentialparam4": 0.0,
        "potentialparam5": 0.0,
        "n_s": {"derived": True},
        "r": {"derived": True},
        "Omega_de": {"derived": True},
        "Omega_m": {"derived": True},
        "sigma8": {"derived": True},
        "age": {"derived": True},
        "rdrag": {"derived": True},
        "zrei": {"derived": True},
        "Y_p": {"derived": True},
        "DH": {"derived": True},
    }

    def initialize(self):
        self.log.info("QICAMBTheory initialized")

    def get_requirements(self):
        """
        Declare sampler-supplied parameters this theory needs.
        """
        return [
            "N_efolds", "gamma", "n_exp",
            "H0", "ombh2", "omch2", "tau", "As", "mnu", "YHe",
            "potentialparam1", "potentialparam2", "potentialparam3", 
            "potentialparam4", "potentialparam5"
        ]

    def calculate(self, state, want_derived=True, **pvs):
        """
        Calculate results for a given set of parameters.
        """
        self.log.debug(f"Entering calculate with want_derived={want_derived}. Initial state keys: {list(state.keys())}")

        camb_params = camb.CAMBparams()
        
        # Inflation parameters
        camb_params.N_efolds = pvs["N_efolds"]
        camb_params.gamma = pvs["gamma"]
        # camb_params.n_exp = pvs["n_exp"] # n_exp is for DE, not inflation
        
        # Cosmological parameters
        camb_params.set_cosmology(
            H0=pvs["H0"],
            ombh2=pvs["ombh2"],
            omch2=pvs["omch2"],
            tau=pvs["tau"],
            mnu=pvs["mnu"],
            YHe=pvs["YHe"],
        )
        
        # Initial power spectrum
        camb_params.InitPower.set_params(As=pvs["As"], ns=0.96) # ns will be calculated

        # Dark energy parameters
        camb_params.DarkEnergy.potentialparam1 = pvs["potentialparam1"]
        camb_params.DarkEnergy.potentialparam2 = pvs["potentialparam2"] # Usually alpha for single exponential
        camb_params.DarkEnergy.potentialparam3 = pvs["n_exp"] # Usually n for single exponential, now mapping to n_exp
        camb_params.DarkEnergy.potentialparam4 = pvs["potentialparam4"] # Usually V1 for double exponential
        camb_params.DarkEnergy.potentialparam5 = pvs["potentialparam5"] # Usually beta for double exponential

        # Specify quintessence model (assuming double exponential here)
        # 0: LambdaCDM, 1: PPF, 2: EarlyDE, 3: Quintessence_std(PL), 4: Quintessence_std(EXP), 
        # 5: Quintessence_std(SIN), 6: Quintessence_std(COSH), 7: Quintessence_std(ASIN), 
        # 8: Quintessence_std(ACOSH), 9: Quintessence_early(EXP)
        camb_params.DarkEnergy.DE_model = 9 # Exponential quintessence

        # We want Cls, and we want lensed Cls
        camb_params.WantCls = True
        camb_params.DoLensing = True
        camb_params.set_for_lmax(lmax=2500, lens_potential_accuracy=1) # lens_potential_accuracy can be tuned
        # Ensure WantTransfer is set for sigma8 calculation
        camb_params.WantTransfer = True

        try:
            self.current_results = camb.get_results(camb_params)
            
            # Store results in state
            state["params"] = pvs
            # Store RAW C_ls directly from CAMB. get_Cl will process these.
            state["_raw_camb_cls"] = self.current_results.get_total_cls(raw_cl=True) 
            state["TCMB"] = camb_params.TCMB # Store TCMB for unit conversion
            
            # Derived parameters
            if want_derived:
                try:
                    derived_params = self.current_results.get_derived_params() # CAMB's own derived params
                    state["n_s"] = self.current_results.Params.InitPower.ns
                    state["r"] = self.current_results.Params.InitPower.r
                    state["Omega_de"] = self.current_results.omega_de # This is Omega_Lambda today
                    # Calculate Omega_m today
                    state["Omega_m"] = self.current_results.get_Omega("cdm", z=0)  + self.current_results.get_Omega("baryon", z=0)
                    state["sigma8"] = self.current_results.get_sigma8_0()
                    state["age"] = derived_params["age"]
                    state["rdrag"] = derived_params["rdrag"]
                    # zrei might not always be present
                    state["zrei"] = derived_params.get("zrei", float('nan')) # Use .get for safety
                    state["Y_p"] = camb_params.YHe # Y_p is approximated by input YHe for now
                    state["DH"] = self.current_results.hubble_parameter(0) / (constants.c / 1000) # H0 in km/s/Mpc, c in m/s. DH wants 1/Mpc.

                except Exception as e:
                    self.log.error(f"Error calculating derived parameters: {e}")
                    import traceback # Make sure traceback is imported
                    self.log.error(traceback.format_exc())
                    # Populate with NaNs for derived params if this block fails
                    state["n_s"] = float('nan')
                    state["r"] = float('nan')
                    state["Omega_de"] = float('nan')
                    state["Omega_m"] = float('nan')
                    state["sigma8"] = float('nan')
                    state["age"] = float('nan')
                    state["rdrag"] = float('nan')
                    state["zrei"] = float('nan')
                    state["Y_p"] = float('nan')
                    state["DH"] = float('nan')
                    # Also ensure _raw_camb_cls and TCMB are present for get_Cl, even if dummy
                    if "_raw_camb_cls" not in state: # Should have been set before this try block
                        dummy_lmax = 2500
                        state["_raw_camb_cls"] = np.zeros((dummy_lmax + 1, 5))
                    if "TCMB" not in state: # Should have been set
                        state["TCMB"] = 2.7255 
                    self.log.error(f"Exiting calculate (derived param error). want_derived={want_derived}. Final state keys: {list(state.keys())}")
                    return False # Signal failure
            
            # If execution reaches here, it means either:
            # 1. want_derived was False (derived params not computed, but that's OK)
            # 2. want_derived was True AND derived params were computed successfully.
            # In both cases, the main CAMB calculation (get_results, get_total_cls) was successful.
            self.log.info(f"Exiting calculate (success). want_derived={want_derived}. Final state keys: {list(state.keys())}")
            return True # Signal success
        
        except camb.CAMBError as e:
            self.log.error(f"CAMB calculation failed: {e}")
            state["logp"] = -np.inf 
            dummy_lmax = 2500
            state["_raw_camb_cls"] = np.zeros((dummy_lmax + 1, 5)) 
            state["TCMB"] = 2.7255 

            if want_derived:
                state["n_s"] = float('nan')
                state["r"] = float('nan')
                state["Omega_de"] = float('nan')
                state["Omega_m"] = float('nan')
                state["sigma8"] = float('nan')
                state["age"] = float('nan')
                state["rdrag"] = float('nan')
                state["zrei"] = float('nan')
                state["Y_p"] = float('nan')
                state["DH"] = float('nan')
            self.log.error(f"Exiting calculate (CAMBError). want_derived={want_derived}. Final state keys: {list(state.keys())}")
            return False 
        
        except Exception as e:
            self.log.error(f"Unexpected error in QICAMBTheory.calculate: {e}")
            import traceback 
            self.log.error(traceback.format_exc())
            state["logp"] = -np.inf
            dummy_lmax = 2500
            state["_raw_camb_cls"] = np.zeros((dummy_lmax + 1, 5))
            state["TCMB"] = 2.7255
            
            if want_derived: # Ensure derived parameters are set to NaN
                state["n_s"] = float('nan')
                state["r"] = float('nan')
                state["Omega_de"] = float('nan')
                state["Omega_m"] = float('nan')
                state["sigma8"] = float('nan')
                state["age"] = float('nan')
                state["rdrag"] = float('nan')
                state["zrei"] = float('nan')
                state["Y_p"] = float('nan')
                state["DH"] = float('nan')
            self.log.error(f"Exiting calculate (unexpected error). want_derived={want_derived}. Final state keys: {list(state.keys())}")
            return False

    def get_can_provide(self):
        """Tell Cobaya what this theory can compute"""
        return ["Cl", "H0", "n_s", "r", "Omega_de", "Omega_m", "sigma8"]

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