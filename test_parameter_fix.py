import camb
from camb import model
import os
import numpy as np

# Create parameters
params = camb.CAMBparams()

# Define output directory
output_dir = 'QI_CAMB/python_test/'
print(f"Attempting to create and use output directory: {os.path.abspath(output_dir)}")
os.makedirs(output_dir, exist_ok=True)

# Use the exact same cosmological parameters from your .ini file
params.set_cosmology(
    H0=67.32117,          # hubble param in km/s/Mpc
    ombh2=0.0223828,      # physical baryon density 
    omch2=0.1201075,      # physical cold dark matter density
    omk=0,                # curvature parameter
    mnu=0.06,             # approximate neutrino mass sum matching omnuh2=0.6451439E-03
    tau=0.05430842        # optical depth
)

# Configure matter power spectrum calculation
params.set_matter_power(redshifts=[0.0], kmax=2.0)

# Define the parameter values to be used
py_V0 = 9e-121
py_alpha = -10.0
py_n = 3.0
py_V1 = 0#8e-123
py_beta = 0# -2.0

phi_evolution_filepath = os.path.join(output_dir, 'phi_evolution_python.dat')
params.set_dark_energy(
    dark_energy_model='early',
    V0=py_V0,
    alpha=py_alpha,
    n=py_n,
    V1=py_V1,
    beta=py_beta,
    output_background_phi=True,
    output_background_phi_filename=phi_evolution_filepath
)

# Print explanation of what we're doing, using the defined variables
print("Setting parameters with DIRECT mapping (values passed to set_dark_energy):")
print(f"  V0 (Python kwarg)    = {py_V0}  -> becomes V0 (param1) in Fortran")
print(f"  alpha (Python kwarg) = {py_alpha}  -> becomes alpha (param2) in Fortran")
print(f"  n (Python kwarg)     = {py_n}       -> becomes n (param3) in Fortran")
print(f"  V1 (Python kwarg)    = {py_V1}      -> becomes V1 (param4) in Fortran")
print(f"  beta (Python kwarg)  = {py_beta} -> becomes beta (param5) in Fortran")

# Set for useful output
params.Want_CMB = True
params.WantTransfer = True
params.DoLensing = True
params.NonLinear = model.NonLinear_both

# Run CAMB
print("\nRunning CAMB...")
try:
    results = camb.get_results(params)
    if results is None:
        print("ERROR: camb.get_results(params) returned None")
    else:
        print("Successfully got results from CAMB.")
        cls = results.get_cmb_power_spectra() 
        print("Success! Model calculated correctly.")

        # Check phi evolution file (saved by Fortran)
        print(f"Checking for phi evolution file at: {os.path.abspath(phi_evolution_filepath)}")
        if os.path.exists(phi_evolution_filepath):
            print(f"Phi evolution file FOUND: {phi_evolution_filepath}")
        else:
            print(f"Phi evolution file NOT FOUND: {phi_evolution_filepath}")

        # Save power spectra results
        cmb_spectra_file = os.path.join(output_dir, 'cmb_power_spectra.dat')
        print(f"Attempting to save CMB power spectra to: {os.path.abspath(cmb_spectra_file)}")
        results.save_cmb_power_spectra(cmb_spectra_file)
        if os.path.exists(cmb_spectra_file):
            print(f"CMB power spectra file FOUND: {cmb_spectra_file}")
        else:
            print(f"CMB power spectra file NOT FOUND: {cmb_spectra_file}")

        #lensing_potential_file = os.path.join(output_dir, 'lensing_potential_cls.dat')
        #print(f"Attempting to save Lensing Potential Cls to: {os.path.abspath(lensing_potential_file)}")
        #results.save_lensing_potential_cls(lensing_potential_file)
        #if os.path.exists(lensing_potential_file):
        #    print(f"Lensing potential Cls file FOUND: {lensing_potential_file}")
        #else:
        #    print(f"Lensing potential Cls file NOT FOUND: {lensing_potential_file}")

        #matter_power_file_path = os.path.join(output_dir, 'matter_power_z0.dat')
        #print(f"Attempting to save Matter Power Spectrum to: {os.path.abspath(matter_power_file_path)}")
        ## Get matter power spectrum data (k, P(k))
        ## For the first redshift params.Transfer.PK_redshifts[0]
        #k, _, PK = results.get_matter_power_spectrum(var1='delta_tot', var2='delta_tot')
        ## Save to file using numpy.savetxt
        #np.savetxt(matter_power_file_path, np.transpose([k, PK]), header='k (Mpc^-1)     P(k) ((Mpc)^3)')
        #if os.path.exists(matter_power_file_path):
        #    print(f"Matter power spectrum file FOUND: {matter_power_file_path}")
        #else:
        #    print(f"Matter power spectrum file NOT FOUND: {matter_power_file_path}")

except Exception as e:
    print("Error:", str(e))