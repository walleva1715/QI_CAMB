from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool, c_char


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int)]

    def validate_params(self):
        return True


class DarkEnergyEqnOfState(DarkEnergyModel):
    """
    Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0.

    See :meth:`.model.CAMBparams.set_initial_power_function` for a convenience constructor function to
    set a general interpolated P(k) model from a python function.

    """
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyEqnOfState'

    _fields_ = [
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("cs2", c_double, "fluid rest-frame sound speed squared"),
        ("use_tabulated_w", c_bool, "using an interpolated tabulated w(a) rather than w, wa above"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)")
    ]

    _methods_ = [('SetWTable', [numpy_1d, numpy_1d, POINTER(c_int)])]

    def set_params(self, w=-1.0, wa=0, cs2=1.0):
        """
         Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed squared
        """
        self.w = w
        self.wa = wa
        self.cs2 = cs2
        self.validate_params()

    def validate_params(self):
        if not self.use_tabulated_w and self.wa + self.w > 0:
            raise CAMBError('dark energy model has w + wa > 0, giving w>0 at high redshift')

    def set_w_a_table(self, a, w):
        """
        Set w(a) from numerical values (used as cublic spline). Note this is quite slow.

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        if len(a) != len(w):
            raise ValueError('Dark energy w(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):
            raise ValueError('Dark energy w(a) arrays must end at a=1')
        if np.any(a <= 0):
            raise ValueError('Dark energy w(a) table cannot be set for a<=0')

        a = np.ascontiguousarray(a, dtype=np.float64)
        w = np.ascontiguousarray(w, dtype=np.float64)

        self.f_SetWTable(a, w, byref(c_int(len(a))))
        return self


@fortran_class
class DarkEnergyFluid(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization using the constant sound-speed single fluid model
    (as for single-field quintessense).

    """

    _fortran_class_module_ = 'DarkEnergyFluid'
    _fortran_class_name_ = 'TDarkEnergyFluid'

    def validate_params(self):
        super().validate_params()
        if not self.use_tabulated_w:
            if self.wa and (self.w < -1 - 1e-6 or 1 + self.w + self.wa < - 1e-6):
                raise CAMBError('fluid dark energy model does not support w crossing -1')


@fortran_class
class DarkEnergyPPF(DarkEnergyEqnOfState):
    """
    Class implementating the w, wa or splined w(a) parameterization in the PPF perturbation approximation
    (`arXiv:0808.3125 <https://arxiv.org/abs/0808.3125>`_)
    Use inherited methods to set parameters or interpolation table.

    """
    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


@fortran_class
class AxionEffectiveFluid(DarkEnergyModel):
    """
    Example implementation of a specifc (early) dark energy fluid model
    (`arXiv:1806.10608 <https://arxiv.org/abs/1806.10608>`_).
    Not well tested, but should serve to demonstrate how to make your own custom classes.
    """
    _fields_ = [
        ("w_n", c_double, "effective equation of state parameter"),
        ("fde_zc", c_double, "energy density fraction at z=zc"),
        ("zc", c_double, "decay transition redshift (not same as peak of energy density fraction)"),
        ("theta_i", c_double, "initial condition field value")]

    _fortran_class_name_ = 'TAxionEffectiveFluid'
    _fortran_class_module_ = 'DarkEnergyFluid'

    def set_params(self, w_n, fde_zc, zc, theta_i=None):
        self.w_n = w_n
        self.fde_zc = fde_zc
        self.zc = zc
        if theta_i is not None:
            self.theta_i = theta_i


# base class for scalar field quintessence models
class Quintessence(DarkEnergyModel):
    r"""
    Abstract base class for single scalar field quintessence models.

    For each model the field value and derivative are stored and splined at sampled scale factor values.

    To implement a new model, need to define a new derived class in Fortran,
    defining Vofphi and setting up initial conditions and interpolation tables (see TEarlyQuintessence as example).

    """
    _fields_ = [
        ("DebugLevel", c_int),
        ("astart", c_double),
        ("integrate_tol", c_double),
        ("sampled_a", AllocatableArrayDouble),
        ("phi_a", AllocatableArrayDouble),
        ("phidot_a", AllocatableArrayDouble),
        ("__npoints_linear", c_int),
        ("__npoints_log", c_int),
        ("__dloga", c_double),
        ("__da", c_double),
        ("__log_astart", c_double),
        ("__max_a_log", c_double),
        ("__ddphi_a", AllocatableArrayDouble),
        ("__ddphidot_a", AllocatableArrayDouble),
        ("omega_tol", c_double),
        ("__state", f_pointer)
    ]
    _fortran_class_module_ = 'Quintessence'


@fortran_class
class EarlyQuintessence(Quintessence):
    r"""
    Double exponential quintessence model with potential

     V(\phi) = V_0 \exp(\alpha \phi^n) + V_1 \exp(\beta \phi)

    This model can accommodate both early dark energy and late-time acceleration.
    """
    _fortran_class_module_ = 'Quintessence'
    _fortran_class_name_ = 'TEarlyQuintessence'
    
    _fields_ = [
        ("n", c_double, "power index for potential"),
        ("f", c_double, r"f/Mpl (sqrt(8\piG)f); legacy parameter"),
        ("m", c_double, "mass parameter; legacy parameter"),
        ("theta_i", c_double, "phi/f initial field value; legacy parameter"),
        ("frac_lambda0", c_double, "fraction of dark energy in cosmological constant today"),
        ("use_zc", c_bool, "legacy parameter"),
        ("zc", c_double, "legacy parameter"),
        ("fde_zc", c_double, "legacy parameter"),
        ("npoints", c_int, "number of points for background integration spacing"),
        ("min_steps_per_osc", c_int, "minimumum number of steps per background oscillation scale"),
        ("fde", AllocatableArrayDouble, "after initialized, the calculated background dark energy "
                                        "fractions at sampled_a"),
        ("__ddfde", AllocatableArrayDouble),
        ("output_background_phi", c_bool, "flag to output background phi evolution"),
        ("output_background_phi_filename", c_char*50, "filename for background phi output"),
        ("search_for_initialphi", c_bool, "flag to search for initial phi"),
        ("potentialparams", c_double * 5, "parameters of the double exponential potential: "
                                           "[V_0, alpha, n, V_1, beta]")
    ]

    def __init__(self, **kwargs):
        # Initialize with default values
        self.search_for_initialphi = False
        
        # Initialize potentialparams to zeros
        for i in range(5):
            self.potentialparams[i] = 0.0
        
        kwargs_copy = kwargs.copy()
        kwargs_copy['search_for_initialphi'] = False
        super().__init__(**kwargs_copy)
        
        self.search_for_initialphi = False

    def set_params(self, V0=1e-120, alpha=0.1, n=1.0, V1=5e-121, beta=0.05,
                 output_background_phi=False, output_background_phi_filename=None):
        """
        Configure double exponential quintessence parameters.
        
        Python keyword arguments map to Fortran potential parameters (P) as follows,
        using a dummy self.potentialparams[0] in Python:
          - Python 'V0' kwarg    -> Fortran P(1) (V0_F)    -> self.potentialparams[0]
          - Python 'alpha' kwarg -> Fortran P(2) (alpha_F) -> self.potentialparams[1]
          - Python 'n' kwarg     -> Fortran P(3) (n_F)     -> self.potentialparams[2]
          - Python 'V1' kwarg    -> Fortran P(4) (V1_F)    -> self.potentialparams[3]
          - Python 'beta' kwarg  -> Fortran P(5) (beta_F)  -> self.potentialparams[4]
        Note: The user's previous setup effectively used self.potentialparams[0] as unused,
        and started mapping from self.potentialparams[1]. This is maintained IF the parameter
        passing in test_parameter_fix.py reflects that desired mapping.
        Given the Fortran output, the mapping seems to be direct from python kwarg name to Fortran param name.
        E.g. Python V0 kwarg sets Fortran V0 (param1).
        This requires Python self.potentialparams[0] = V0_kwarg, etc.
        """
        self.search_for_initialphi = False
        
        for i in range(5):
            self.potentialparams[i] = 0.0
            
        # Direct mapping: Python kwarg name corresponds to Fortran parameter name.
        # Python self.potentialparams is 0-indexed, Fortran potentialparams is 1-indexed.
        # So, self.potentialparams[0] maps to Fortran potentialparams(1) (V0_F)
        # self.potentialparams[1] maps to Fortran potentialparams(2) (alpha_F)
        # ...and so on.
        self.potentialparams[0] = V0     # Python V0 kwarg -> Fortran P(1) (V0 in formula)
        self.potentialparams[1] = alpha  # Python alpha kwarg -> Fortran P(2) (alpha in formula)
        self.potentialparams[2] = n      # Python n kwarg -> Fortran P(3) (n in formula)
        self.potentialparams[3] = V1     # Python V1 kwarg -> Fortran P(4) (V1 in formula)
        self.potentialparams[4] = beta   # Python beta kwarg -> Fortran P(5) (beta in formula)
        
        self.output_background_phi = output_background_phi
        if output_background_phi_filename is not None:
            # Ensure string fits in the 50-char buffer and convert to bytes
            fname = str(output_background_phi_filename)[:50]
            self.output_background_phi_filename = fname.encode('utf-8')
        
        self.search_for_initialphi = False
        
        return self


# short names for models that support w/wa
F2003Class._class_names.update({
    'fluid': DarkEnergyFluid,
    'ppf': DarkEnergyPPF,
    'early': EarlyQuintessence
})
