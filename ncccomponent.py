''' OpenMDAO component wrapper for NCC '''
# --- Inherent python/system level imports
import os
import math
import sys
import logging 
import shutil
import datetime
from glob import iglob

# --- OpenMDAO main and library imports
from openmdao.main.api import Component, FileMetadata, enable_console
from openmdao.lib.datatypes.api import Enum, Float, Int, Str, Bool, Event
from openmdao.lib.components.api import ExternalCode
from openmdao.main.resource import ResourceAllocationManager as RAM
from nas_access import NAS_Allocator

# --- OpenMDAO component imports
from pathname import Path
        
def _bool2str(var):
    ''' Converts OpenMDAO boolean variable to string,
    'T' if True, 'F' if False. '''
    
    if var == True:
        return 'T'.center(25)
    else:
        return 'F'.center(25)
        
def _int2str(var):
    ''' Converts OpenMDAO integer variable to string '''
    
    return str(var).center(25)
    
def _enum2str(var):
    ''' Converts OpenMDAO enumerated variable to string '''
    
    return str(var).center(25) 

def copy_files(src_glob, dst_folder):
    for fname in iglob(src_glob):
        shutil.copy2(fname, dst_folder)    
    
# ==================================================================================================================
# ===                           File Wrapper for the National Combustion Code                                    ===
# ==================================================================================================================
class NCCcomponent(ExternalCode):
    ''' OpenMDAO component wrapper for NCC '''
    
    clean_start = Event()
    
    def _clean_start_fired(self):
        self._clean_start = True
    
# ==================================================================================================================
# ===                                       Initialize Defaults                                                  ===
# ==================================================================================================================        
    
    config = Str('1', iotype = 'in', desc = 'current design configuration number')     
    low_dissipation = Bool(False, iotype = 'in', desc = 'specify whether dissipation level is low')
    FAR = Float(0.0, iotype = 'in', desc = 'total combustor fuel-to-air ratio')
    rho_fuel = Float(763.0, iotype = 'in', desc = 'fuel density', units = 'kg/m**3')
    rho_air = Float(14.2649, iotype = 'in', desc = 'combustor air density', units = 'kg/m**3')  
    
    num_injectors = Int(5, iotype = 'in', desc = 'number of fuel injector modules that are active')
    injector_dia = Float(0.0, iotype = 'in', desc = 'injector diameter')
    
    nonreacting = Bool(False, iotype = 'in', desc = 'flag specifier to run ncc nonreacting')
    reacting = Bool(False, iotype = 'in', desc = 'flag specifier to run ncc with liquid spray and artificial ignition')
    reacting2 = Bool(False, iotype = 'in', desc = 'flag specifier to run reacting with detailed mechanism')    
    
    # ------------------------------------------
    # --- Initialize values for ncc_input.in ---
    # ------------------------------------------
    axisymmetric = Bool(False, iotype = 'in', desc = 'case axisymmetric T/F')
    compressible = Bool(False, iotype = 'in', desc = 'case compressible T/F')
    euler = Bool(False, iotype = 'in', desc = 'apply inviscid euler equations T/F')
    time_accurate = Bool(False, iotype = 'in', desc = 'LES T/F')
    real_time_step = Float(0.0, iotype = 'in', desc = 'step size')
    
    continuity_goal = Int(2000, iotype = 'in', desc = '# of consecutive iterations in which mass imbalence goal is met')
    mass_imbalance_goal = Float(1.0e-3, iotype = 'in', desc = '(m_in - m_out)/m_in')
    CFL = Float(1.0, iotype = 'in', desc = 'Dictates rate of convergence by adjusting step size')
    max_multigrid_level = Int(1, iotype = 'in')
    
    max_iterations_per_time_step = Int(50000, iotype = 'in', desc = 'max number of iterations')
    real_time_iterations = Int(1, iotype = 'in', desc = 'number of real time iterations')
    restarts = Int(50, iotype = 'in', desc = 'max number of restarts')
    
    turbulence = Bool(True, iotype = 'in', desc = 'model turbulence T/F')
    low_re_model = Bool(False, iotype = 'in', desc = 'apply low Re number model')
    var_cmu = Bool(True, iotype = 'in')
    quadratic = Bool(True, iotype = 'in')
    cubic = Bool(True, iotype = 'in')
    scalar_nl = Bool(False, iotype = 'in')
    spec_nl = Bool(False, iotype = 'in')
    species = Bool(True, iotype = 'in', desc = 'model species')
    enthalpy = Bool(True, iotype = 'in', desc = 'calculate enthalpy')
    
    mean_on = Bool(False, iotype = 'in', desc = 'turn on for unsteady simulations')
    residual_smoothing_on = Bool(True, iotype = 'in', desc = 'turn on for steady simulations')
    
    iter_turbulence = Int(1, iotype = 'in', desc = 'turbulence iterations')
    iter_species = Int(1, iotype = 'in', desc = 'species iterations')
    iter_enthalpy = Int(1, iotype = 'in', desc = 'enthalpy iterations')

    inv_sigma_k = Float(1.00, iotype = 'in', desc = '1/sigma_k')
    inv_sigma_e = Float(0.769231, iotype = 'in', desc = '1/sigma_e')
    
    RCp = Float(1.0, iotype = 'in', desc = 'use < 1 for TFNS')
    C_mu = Float(0.09, iotype = 'in')
    C_e1 = Float(1.45, iotype = 'in')
    C_e2 = Float(1.92, iotype = 'in')
    
    laminar_viscosity = Float(1.835e-5, iotype = 'in')
    effect_lam_visc = Float(1.0e4, iotype = 'in')
    under_relax = Float(0.9, iotype = 'in')
    RNG_cap = Float(5.5, iotype = 'in')
    
    inv_Sch_1 = Float(1.43, iotype = 'in')
    inv_Sch_2 = Float(1.43, iotype = 'in')
    inv_Sch_3 = Float(1.43, iotype = 'in')
    inv_Sch_4 = Float(1.43, iotype = 'in')
    
    stoic = Float(1.0, iotype = 'in')
    A1_EDC = Float(4.0, iotype = 'in')
    A2_EDC = Float(2.0, iotype = 'in')
  
    universal_gas_constant = Float(8314.0, iotype = 'in')
    mw_fuel = Float(28.97, iotype = 'in')
    mw_oxidant = Float(28.97, iotype = 'in')
    mw_prod = Float(28.97, iotype = 'in')
 
    cp_fuel = Float(1005.0, iotype = 'in')
    cp_oxidant = Float(1005.0, iotype = 'in')
    cp_prod = Float(1005.0, iotype = 'in')
    Prandtl = Float(0.713, iotype = 'in')
 
    combust = Bool(False, iotype = 'in') 
    chem_model = Enum(1, [0, 1, 2], iotype = 'in')
    href_fuel = Float(0.0, iotype = 'in')
    href_oxidant = Float(0.0, iotype = 'in')
    href_prod = Float(0.0, iotype = 'in')
    Tref = Float(298.15, iotype = 'in')
   
    min_temperature = Float(250.0, iotype = 'in')
    max_temperature = Float(3200.0, iotype = 'in')
    
    lspray = Bool(False, iotype = 'in')
    num_liquids = Int(0, iotype = 'in', desc = 'number of fuel injectors')
    
    aero2 = Float(-1E-3, iotype = 'in')
    k_e2 = Float(-1E-3, iotype = 'in')
    species2 = Float(-1E-3, iotype = 'in')
    enthalpy2 = Float(-1E-3, iotype = 'in')
    
    aero4 = Float(0.05, iotype = 'in')
    k_e4 = Float(0.05, iotype = 'in')
    species4 = Float(0.05, iotype = 'in')
    enthalpy4 = Float(0.05, iotype = 'in')

    aero4_switch = Bool(True, iotype = 'in')
    k_e4_switch = Bool(True, iotype = 'in')
    species4_switch = Bool(True, iotype = 'in')
    enthalpy4_switch = Bool(True, iotype = 'in')
    
    k2 = Float(0.5, iotype = 'in')                  
    beta_scale = Float(1.0, iotype = 'in')
    
    expn = Float(0.7, iotype = 'in')
    
    high_parallel = Bool(True, iotype = 'in')

    chem_turb_model = Int(0, iotype = 'in') 
    rxn_low_temp = Float(500.0, iotype = 'in')
    
    y_fuel_fuel = Float(0.0, iotype = 'in')
    y_o2_fuel = Float(0.0, iotype = 'in')
    y_co2_fuel = Float(0.0, iotype = 'in')
    y_h2o_fuel = Float(0.0, iotype = 'in')
    y_n2_fuel = Float(0.0, iotype = 'in')
    
    y_fuel_oxid = Float(0.0, iotype = 'in')
    y_o2_oxid = Float(0.0, iotype = 'in')
    y_co2_oxid = Float(0.0, iotype = 'in')
    y_h2o_oxid = Float(0.0, iotype = 'in')
    y_n2_oxid = Float(0.0, iotype = 'in')
    
    debug_res_on = Bool(True, iotype = 'in')
    debug_res_start = Int(0, iotype = 'in') 
    debug_res_print = Int(10, iotype = 'in') 

    pg_on = Bool(False, iotype = 'in')

    hp_eos = Int(0, iotype = 'in') 
    
    vnn_solid = Float(2.0, iotype = 'in')

    angle_axi = Float(2.0, iotype = 'in')

    int_chem_solver = Int(-1, iotype = 'in')
    substeps = Int(10, iotype = 'in')
    
    newton_temp = Float(6.0, iotype = 'in')
    tol_temp = Float(0.01, iotype = 'in')

    calc_lam_diff = Bool(False, iotype = 'in')
    mass_diff_in_heat_flux = Bool(False, iotype = 'in')
    
    max_iter_tstatic = Float(6.0, iotype = 'in')
    newton_tstatic = Float(6.0, iotype = 'in')
    tol_tstatic = Float(0.05, iotype = 'in')

    over_write = Bool(True, iotype = 'in')
    add_mut_switch = Bool(False, iotype = 'in')
    
    ce3 = Float(1.0, iotype = 'in')
    add_ce3_term = Bool(False, iotype = 'in')

    enth_variance_on = Bool(False, iotype = 'in')
    iter_enth_variance = Int(1, iotype = 'in')
    
    tprime_intensity = Float(0.001, iotype = 'in')
    g_decrease_factor = Float(0.9, iotype = 'in')

    apdf_tprime_small = Float(10.0, iotype = 'in')
    apdf_tprime_large = Float(1000.0, iotype = 'in')

    apdf_simpson_eps = Float(1.0e-03, iotype = 'in')
    apdf_simpson_jmax = Int(10, iotype = 'in')

    thermal_wallfct = Bool(False, iotype = 'in')

    WRITE_DB_CONTROL = Int(1, iotype = 'in')

    # --------------------------------------------
    # --- Initialize values for ncc_bcflags.in ---
    # --------------------------------------------
    bc1_mdot = Float(0.0, iotype = 'in', desc = 'inlet air velocity', units = 'kg/s')
    bc1_Tstatic = Float(700, iotype = 'in', desc = 'inlet air static temperature', units = 'degK')    
    bc1_turb = Float(0.05, iotype = 'in', desc = 'inlet air turbulence')
    bc1_Lmix = Float(0.0254, iotype = 'in', desc = 'inlet air mixing length (typically 1/3 flowpath diameter) ', units = 'm')
    bc1_Mach_no = Float(0.04, iotype = 'in', desc = 'inlet air Mach number guess')    
    bc1_O2 = Float(0.23, iotype = 'in')
    bc1_N2 = Float(0.77, iotype = 'in')   
    
    bc2_Pstatic = Float(101325.0, iotype = 'in', desc = 'exit static pressure', units = 'Pa')    
    bc2_Tstatic = Float(700.0, iotype = 'in', desc = 'exit static temperature', units = 'degK')
    bc2_subsonic = Bool(True, iotype = 'in', desc = 'exit subsonic T/F')

    bc7_model = Int(2, iotype = 'in', desc = 'model type (2 = no-slip)')
    bc7_Twall = Float(-300.0, iotype = 'in', desc = 'wall temperature (< 0 = adiabatic)', units = 'degK')
    bc7_mdot = Float(0.0, iotype = 'in', desc = 'cooling air velocity', units = 'kg/s')
    bc7_Aflow = Float(0.000388, iotype = 'in', desc = 'area of effusion holes', units = 'm**2')
    bc7_Tstatic = Float(700.0, iotype = 'in', desc = 'cooling air velocity', units = 'degK')
    bc7_turb = Float(0.05, iotype = 'in', desc = 'cooling air turbulence')
    bc7_Lmix = Float(0.00125, iotype = 'in', desc = 'cooling air mixing length (typically 1/3 flowpath diameter)', units = 'm')   
 
    # ---------------------------------------------------
    # --- Initialize values for ncc_material_flags.in ---
    # ---------------------------------------------------
    axisymmetric_rotation = Float(0.0, iotype = 'in')
    
    u_init = Float(0.0, iotype = 'in', desc = 'initial axial velocity', units = 'm/s')
    v_init = Float(0.0, iotype = 'in', desc = 'initial horizontal velocity', units = 'm/s')
    w_init = Float(0.0, iotype = 'in', desc = 'initial vertical velocity', units = 'm/s')
    Pstatic_init = Float(0.0, iotype = 'in', desc = 'initial static pressure', units = 'Pa')
    k_init = Float(0.0, iotype = 'in')
    e_init = Float(0.0, iotype = 'in')
    Tstatic_init = Float(700.0, iotype = 'in', desc = 'initial static temperature', units = 'degK')
    
    O2_init = Float(0.23, iotype = 'in')
    N2_init = Float(0.77, iotype = 'in')
    
    # ----------------------------------------------
    # --- Initialize values for nodal_results.in ---
    # --- Used to calculate variables for Tecplot --
    # ----------------------------------------------
    file_num = Int(1, iotype = 'in')
    results_fn = Str('results', iotype = 'in')              
        
    # --------------------------------------------------
    # --- Initialize values for ncc_walltime_stop.in ---
    # --------------------------------------------------
    twall_run = Float(7.0, iotype = 'in')
    twall_restart_files = Float(30.0, iotype = 'in')
    twall_factor = Float(1.0, iotype = 'in') 
  
    # ---------------------------------------------
    # --- Initialize values for ncc_ignition.in ---
    # ---------------------------------------------    
    ignition_on = Bool(False, iotype = 'in')
    iter_ignition = Int(-1, iotype = 'in')
    no_of_ignite = Int(1, iotype = 'in')
   
    xmin1 = Float(0.0, iotype = 'in')
    xmax1 = Float(0.0, iotype = 'in')
    ymin1 = Float(0.0, iotype = 'in')
    ymax1 = Float(0.0, iotype = 'in')
    zmin1 = Float(0.0, iotype = 'in')
    zmax1 = Float(0.0, iotype = 'in')

    limit_ignition_duration = Int(3000, iotype = 'in')
    
    product_seed_value = Float(1.0e-2, iotype = 'in')
    prod_off = Float(1.1, iotype = 'in')
    
    delT = Int(5, iotype = 'in')
    T_off = Int(1800, iotype = 'in')
    ignite_on_restart = Int(-1, iotype = 'in')
    spark_cap = Float(1.0e12, iotype = 'in')
    
    ignite_done_model = Int(1, iotype = 'in')
    ignite_off_model = Int(1, iotype = 'in')
    inner_radius = Float(-0.007, iotype = 'in')

    ncells_off = Int(10, iotype = 'in')
    vol_off = Float(0.0, iotype = 'in')
    
    energy_ignite_total = Float(0.0, iotype = 'in')
    time_ignite_total = Float(0.0, iotype = 'in')    
    
    spec_name_ignite = Str('C11H21', iotype = 'in')    
    
    del_yisp = Float(0.0, iotype = 'in')
    yisp_off = Float(0.0, iotype = 'in')  
    fact_ILDM_ig = Float(0.1, iotype = 'in')
    
    # -----------------------------------------------
    # --- Initialize values for ncc_injector.in.# ---
    # -----------------------------------------------    
    number_of_liquid_components = Int(1, iotype = 'in')
    
    fuel_symbol = Str('C11H21', iotype = 'in')
    
    mass_fraction = Float(1.0, iotype = 'in')
    
    initial_liquid_temperature = Float(315.0, iotype = 'in')   
    superheat_vap_model = Int(0, iotype = 'in')
    Jet_A_properties = Bool(False, iotype = 'in')
    init_mass_fraction_vapor = Float(0.0, iotype = 'in')
    
    atomization = Bool(False, iotype = 'in')
    drop_breakup_model = Bool(False, iotype = 'in')
    line_injection = Bool(False, iotype = 'in')
    spray_table = Bool(False, iotype = 'in')
    steady_spray_model = Bool(False, iotype = 'in')    
    
    no_of_holes = Int(1, iotype = 'in') 
    no_of_streams = Int(16, iotype = 'in')
    no_of_droplet_groups = Int(10, iotype = 'in')
    
    lmdis = Int(1, iotype = 'in') 
    smdm = Float(14.0, iotype = 'in')    
    cone = Bool(True, iotype = 'in')
    size_min = Float(1.0, iotype = 'in')
    size_max = Float(100.0, iotype = 'in')
    stochastic_injection = Bool(True, iotype = 'in')
    
    x = Float(0.0, iotype = 'in')
    y = Float(0.0, iotype = 'in')
    z = Float(0.0, iotype = 'in')
    mdot_spray = Float(0.0, iotype = 'in')
    v_inj = Float(25.0, iotype = 'in')
    alpha = Float(0.0, iotype = 'in')
    beta = Float(0.0, iotype = 'in')
    theta = Float(70.0, iotype = 'in')
    dtheta = Float(15.0, iotype = 'in')
    swrl_angle = Float(0.0, iotype = 'in')
    
    atom_type = Int(0, iotype = 'in')
    breakup_type = Int(0, iotype = 'in')
    dia_hole = Float(1.0e-4, iotype = 'in')
    delp_inj = Float(400.0, iotype = 'in')
    liq_vel = Float(25.0, iotype = 'in')
    gas_u = Float(0.0, iotype = 'in')
    gas_v = Float(0.0, iotype = 'in')
    gas_w = Float(0.0, iotype = 'in')
    pcl_start = Int(1, iotype = 'in')
    pcl_end = Int(24, iotype = 'in')
    
    line_injection_type = Int(1, iotype = 'in')
    xi = Float(1.0, iotype = 'in')
    yi = Float(1.0, iotype = 'in')
    zi = Float(1.0, iotype = 'in')
    xf = Float(1.0, iotype = 'in')
    yf = Float(1.0, iotype = 'in')
    zf = Float(1.0, iotype = 'in')
    ri = Float(1.0, iotype = 'in')
    
    # --------------------------------------------------
    # --- Initialize values for ncc_liquid_solver.in ---
    # --------------------------------------------------
    ldread = Bool(True, iotype = 'in')
    ispray_mod = Int(5, iotype = 'in')
    when_start_spray = Int(1, iotype = 'in')

    iswim = Enum(0, [0, 1], iotype = 'in')
    wall_film_thickness = Float(0.0, iotype = 'in')
    offset = Float(0.0, iotype = 'in')
    coeff_of_restitution = Float(2.0, iotype = 'in')
    
    dtml = Float(0.2e-6, iotype = 'in')
    dtgl = Float(1.0e-6, iotype = 'in')
    dtil = Float(1.0e-6, iotype = 'in')
    
    evap_calc_mod = Int(1, iotype = 'in')
    
    sklim_constant = Float(0.04, iotype = 'in')
    
    use_constant_rholm = Bool(True, iotype = 'in')
    
    tol_bmgs = Float(1.0e-6, iotype = 'in')
    
    tol_xtold_lower = Float(0.5, iotype = 'in')
    tol_xtold_upper = Float(1.5, iotype = 'in')
    
    tol_ymk_lower = Float(0.75, iotype = 'in')
    tol_ymk_upper = Float(1.25, iotype = 'in')
    
    tol_xt_tdrop = Float(0.75, iotype = 'in')
    tol_xt_tboil = Float(0.9975, iotype = 'in')
    
    force_evap_mass_factor = Float(1.0, iotype = 'in')
    
    num_calls_no_evap = Int(1, iotype = 'in')
    
    evaporation_on = Bool(True, iotype = 'in')
    
    i_interp_gas_part = Int(1, iotype = 'in')
    
    tol_prob_bprime_sprips = Float(-0.99, iotype = 'in')
    bprime_chas_1 = Float(-0.99, iotype = 'in')
    bprime_chas_2 = Float(-0.99, iotype = 'in')

    tol_prob_devu1_sprips = Float(0.00005, iotype = 'in')
    devu1_chas_1 = Float(0.000001, iotype = 'in')
    devu1_chas_2 = Float(0.000001, iotype = 'in')
    
    tol_prob_sklim_1 = Float(0.1, iotype = 'in')
    sklim_2 = Float(0.1, iotype = 'in')
    
    tol_prob_speed_chg_1 = Float(1.25, iotype = 'in')
    chg_2 = Float(2.0, iotype = 'in')
    
    relax_gas_src_mass = Float(1.0, iotype = 'in')
    moment = Float(1.0, iotype = 'in')
    energy = Float(1.0, iotype = 'in')
    
    i_relax_gas_src_spray = Int(1, iotype = 'in')
    
    max_spray_int_count = Int(9000, iotype = 'in')
    spray_int_interval = Int(5, iotype = 'in')
    
    i_update_icrit = Int(1, iotype = 'in')
    
    stochastic_spray = Bool(True, iotype = 'in')
    
    max_part_searches = Int(-99, iotype = 'in')
    max_part_search_exchanges = Int(2000, iotype = 'in')
    
    limit_dtml = Bool(False, iotype = 'in')                                                             
    ncells_per_dtml = Float(1.0, iotype = 'in')
    
    keyword_partsearch = Str('macp', iotype = 'in')
    rklim = Float(1E-10, iotype = 'in')
       
    # --------------------------------------------------
    # --- Initialize values for ncc_pm_input.in ---
    # --------------------------------------------------
    pbcfl = Float(1.0, iotype = 'in', desc = 'local CFL gradient limiter')
    etau_beta = Float(0.15, iotype = 'in', desc = 'local energy gradient limiter')

# ================================================================================================    
    def __init__(self, *args, **kwargs):
        '''
         Constructor for the NCC component ---
        '''
        super(NCCcomponent, self).__init__(*args, **kwargs)         
        
        self._clean_start = False
        self._config = ''
        self._num_procs = 400
        self.force_execute = True 
        
# ==================================================================================================================
# ===                                               Execution Block                                              ===
# ==================================================================================================================
    def execute(self):
        '''
         Executes file-wrapped component --- 
        '''          

        copy_files(os.path.join(Path('Scripts'), '*'), Path('OpenMDAO'))

        walltime = str(datetime.timedelta(seconds=self.twall_run*3600))

        queue = 'normal' 

        if (self.twall_run <= 2.0):  # --- Allow for faster turn-around time on short jobs
            queue = 'devel'

        if self.nonreacting:
            self._config = self.config + '_NR'                
            self.command = ['python', 'qsub.py', '/run_nonreacting', str(self._num_procs), walltime, queue, 'Sim-' + self._config]
        elif self.reacting:
            self._config = self.config + '_R'
            self.command = ['python', 'qsub.py', '/run_reacting', str(self._num_procs), walltime, queue, 'Sim-' + self._config]                                    
        else:          
            self._config = self.config + '_R'
            self.command = ['python', 'qsub.py', '/run_reacting2', str(self._num_procs), walltime, queue, 'Sim-' + self._config]

            # ---------------------------------------
            # --- Read Previous Execution Results ---
            # ---------------------------------------
            prev_config = self._config.split('_')
            prev_config[1] = str(int(prev_config[1]) - 1)

            prev_dir = os.path.join(Path('Tecplot'), 'Sim-' +  '_'.join(prev_config), 'Input')

            try:
                with open(os.path.join(prev_dir, 'ncc_pm_input.in'), 'r') as f:
                    self.etau_beta = float(f.readline().split()[-1])
            except:
                pass

            self.pbcfl = 0.5
            self.etau_beta = self.etau_beta - 0.05
            self.energy = self.energy*10.0

            if self.energy > 1.0: self.energy = 1.0
            if self.etau_beta < 0.001: self.etau_beta = 0.001

        # ------------------------------------------------------------------------
        # --- Initialize values for ncc_chemistry.in 1-Step Reaction Mechanism ---
        # ------------------------------------------------------------------------
        if self.nonreacting:
            line1 = 'ELEMENTS C H O N END'
            line2 = 'SPECIES C12H23 O2 CO2 H2O N2 END'
            line3 = 'FUEL C12H23 END'
            line4 = 'OXID O2     END'
            line5 = 'REACTIONS'
            line6 = '4 C12H23  + 71 O2  => 48 CO2  + 46 H2O           8.60E+11  0.00  3.00E+4'
            line7 = 'GLO / C12H23 0.10 /'
            line8 = 'GLO / O2     1.65 /'
            line9 = 'END' 
                       
            # -----------------------------------------
            # --- Construct ncc_chemistry data list ---
            # ----------------------------------------- 
            ncc_chemistry_data = [line1]
            ncc_chemistry_data.append(line2)            
            ncc_chemistry_data.append(line3)
            ncc_chemistry_data.append(line4)
            ncc_chemistry_data.append(line5)
            ncc_chemistry_data.append(line6)
            ncc_chemistry_data.append(line7)
            ncc_chemistry_data.append(line8)
            ncc_chemistry_data.append(line9)                        
        
        else:                  
            
            if self.reacting:
                # -------------------------------------------------------------------------
                # --- Initialize values for ncc_chemistry.in 1-Step Reaction Mechanism ---
                # -------------------------------------------------------------------------            
                line1 = 'ELEMENTS C H O N END'
                line2 = 'SPECIES C12H23 O2 CO2 H2O N2 END'
                line3 = 'FUEL C12H23 END'
                line4 = 'OXID O2     END'
                line5 = 'REACTIONS'
                line6 = '4 C12H23  + 71 O2  => 48 CO2  + 46 H2O           8.60E+11  0.00  3.00E+4'
                line7 = 'GLO / C12H23 0.10 /'
                line8 = 'GLO / O2     1.65 /'
                line9 = 'END' 
                
                # -----------------------------------------
                # --- Construct ncc_chemistry data list ---
                # ----------------------------------------- 
                ncc_chemistry_data = [line1]
                ncc_chemistry_data.append(line2)            
                ncc_chemistry_data.append(line3)
                ncc_chemistry_data.append(line4)
                ncc_chemistry_data.append(line5)
                ncc_chemistry_data.append(line6)
                ncc_chemistry_data.append(line7)
                ncc_chemistry_data.append(line8)
                ncc_chemistry_data.append(line9)            
        
            else:             
                # -------------------------------------------------------------------------
                # --- Initialize values for ncc_chemistry.in 18-Step Reaction Mechanism ---
                # ------------------------------------------------------------------------- 
                line1 = 'ELEMENTS  C H O N  END'
                line2 = 'SPECIES C11H21 CO2 CO O2 O OH N2 H2 H H2O HO2 NO N2O N CH END'
                line3 = 'REACTIONS'
                line4 = 'C11H21   +   O2  => 11CH  + 10H + O2      1.00E+12  0.       ' + str(round(30192.0 * math.log(self.bc1_Tstatic) - 169422.0, 4))  
                line5 = 'GLO / C11H21 0.8 /'
                line6 = 'GLO / O2     0.8 /'
                line7 = 'CH   + O2        => CO + OH               2.00E+15  0.    0.'
                line8 = 'CH   + O         => CO + H                3.00E+12  1.00  0.'            
                line9 = 'H2   + O2       <=> H2O  + O              3.98E+11  1.00  4.80E+4'
                line10 = 'H2   + O        <=> H    + OH             3.00E+14  0.    6.00E+3'
                line11 = 'H    + O2       <=> O    + OH             4.00E+14  0.    1.80E+4'
                line12 = 'CO   + OH       <=> CO2  + H              5.51E+07  1.27 -7.58E+2'
                line13 = 'H2O  + O2       <=> 2O   + H2O            3.17E+12  2.00  1.12E+5'
                line14 = 'CO   + H2O      <=> CO2  + H2             5.50E+04  1.28 -1.00E+3'
                line15 = 'CO   + H2+O2    <=> CO2  + H2O            1.60E+14  1.60  1.80E+4'
                line16 = 'N    + NO       <=> N2   + O              3.00E+12  0.30  0.'
                line17 = 'N    + O2       <=> NO   + O              6.40E+09  1.00  3.17221E+3'
                line18 = 'N    + OH       <=> NO   + H              6.30E+11  0.50  0.'   
                line19 = 'N    + N + M    <=> N2   + M              2.80E+17 -0.75  0.'
                line20 = 'H    + N2O      <=> N2   + OH             3.50E+14  0.    755.29'
                line21 = 'N2   + O2 + O   <=> N2O  + O2             1.00E+15  0.    302.1'
                line22 = 'N2O  + O        <=> 2NO                   1.50E+15  0.    39000.'
                line23 = 'N2O  + M        <=> N2   + O + M          1.16E+15  0.    33232.63'            
                line24 = 'END'              
                  
                # -----------------------------------------
                # --- Construct ncc_chemistry data list ---
                # ----------------------------------------- 
                ncc_chemistry_data = [line1]
                ncc_chemistry_data.append(line2)            
                ncc_chemistry_data.append(line3)
                ncc_chemistry_data.append(line4)
                ncc_chemistry_data.append(line5)
                ncc_chemistry_data.append(line6)
                ncc_chemistry_data.append(line7)
                ncc_chemistry_data.append(line8)
                ncc_chemistry_data.append(line9) 
                ncc_chemistry_data.append(line10)
                ncc_chemistry_data.append(line11)
                ncc_chemistry_data.append(line12)
                ncc_chemistry_data.append(line13)
                ncc_chemistry_data.append(line14)
                ncc_chemistry_data.append(line15)
                ncc_chemistry_data.append(line16)
                ncc_chemistry_data.append(line17)
                ncc_chemistry_data.append(line18)
                ncc_chemistry_data.append(line19)
                ncc_chemistry_data.append(line20)
                ncc_chemistry_data.append(line21)
                ncc_chemistry_data.append(line22)
                ncc_chemistry_data.append(line23)
                ncc_chemistry_data.append(line24)                         
      
            self.num_liquids = self.num_injectors
            
            # -----------------------------------------------        
            # --- Calculate liquid spray model parameters ---
            # -----------------------------------------------            
            self.mdot_spray = self.FAR * self.bc1_mdot / self.num_injectors          
            kinematic_viscosity = 0.87E-6 # m^2/s @ 80 C
            dynamic_viscosity = kinematic_viscosity * self.rho_fuel # --- kg /(m*s)
            surface_tension = 0.0186 # --- N/m @ 80 C
            percent_air_core = 0.2
            fuel_area = self.mdot_spray / (self.v_inj * self.rho_fuel) # --- m^2
            air_core_area = percent_air_core * fuel_area # --- m^2
            orifice_area = fuel_area + air_core_area
            dP_fuel = 100 # --- psi
                        
            self.smdm = 2.25 * surface_tension**0.25 * dynamic_viscosity**0.25 * self.mdot_spray**0.25 * \
                        (dP_fuel * 6894.75729)**-0.5 * self.rho_air**-0.25 * 1E6
            
            self.v_inj = (10403 * self.mdot_spray - 0.3151)
            self.liq_vel = self.v_inj
            self.smdm = 8.8
            
            print '------------------------------------------------'
            print '---          Injector Parameters             ---'
            print '------------------------------------------------'            
            print 'number injectors    = ', self.num_injectors
            print 'FAR                 = ', self.FAR        
            print 'mdot air            = ', self.bc1_mdot, ' kg/s'
            print 'mdot spray          = ', self.mdot_spray, 'kg/s'        
            print 'fuel density        = ', self.rho_fuel, ' kg/m3' # --- @ 325 K (51.85 C)
            print 'dynamic viscosity   = ', dynamic_viscosity, ' kg/ms' # --- @ 325 K (51.85 C)
            print 'kinematic viscosity = ', kinematic_viscosity, ' m^2/s' # --- @ 325 K (51.85 C)        
            print 'air density         = ', self.rho_air, 'kg/m^3'                
            print 'spray velocity      = ', self.v_inj, ' m/s' 
            print 'SMD                 = ', self.smdm, ' microns'
            print '-------------------------------------------------'                    

        # ------------------------------------------
        # --- Get injector and ignition location ---
        # ------------------------------------------                              
        injector_loc = []
        ignition_loc = []

        line_num = 0

        with open(Path('OpenMDAO') + '\inj_loc.txt', 'r') as f:

            while line_num < self.num_injectors*2:
                line = f.readline()
                
                if line_num % 2 == 0:
                    injector_loc.append(line.split())
                else:
                    ignition_loc.append(line.split())
                line_num+=1
        
        # ==============================================================================================================
        # ===                               Construct Input file data lists                                          ===
        # ==============================================================================================================            
        # -------------------------------------
        # --- Construct ncc_input data list --- 
        # -------------------------------------        
        ncc_input_data = ['Title : Combustor Design Simulation']
        ncc_input_data.append(' '.join(['axisymmetric'.center(25), 'compressible'.center(25), 'euler'.center(25), 'time_accurate'.center(25), 'real_time_step'.center(25)]))
        ncc_input_data.append('{0:25s} {1:25s} {2:25s} {3:25s} {4:25.16f}'.format(_bool2str(self.axisymmetric), _bool2str(self.compressible), _bool2str(self.euler), _bool2str(self.time_accurate), self.real_time_step))
        ncc_input_data.append(' '.join(['continuity_goal'.center(25), 'mass_imbalance_goal'.center(25), 'CFL'.center(25), 'max_multigrid_level'.center(25)]))            
        ncc_input_data.append('{0:25s} {1:25.16f} {2:25.16f} {3:25s}'.format(_int2str(self.continuity_goal), self.mass_imbalance_goal, self.CFL, _int2str(self.max_multigrid_level)))        
        ncc_input_data.append(' '.join(['max_iter_per_tstep'.center(25), 'real_time_iterations'.center(25), 'restarts'.center(25)]))        
        ncc_input_data.append('{0:25s} {1:25s} {2:25s}'.format(_int2str(self.max_iterations_per_time_step), _int2str(self.real_time_iterations), _int2str(self.restarts)))         
        ncc_input_data.append(' '.join(['turbulence'.center(25), 'low_re_model'.center(25), 'var_cmu'.center(25), 'quadratic'.center(25), 'cubic'.center(25), 'scalar_nl'.center(25), 'spec_nl'.center(25), 'species'.center(25), 'enthalpy'.center(25)]))
        ncc_input_data.append('{0:25s} {1:25s} {2:25s} {3:25s} {4:25s} {5:25s} {6:25s} {7:25s} {8:25s}'.format(_bool2str(self.turbulence), _bool2str(self.low_re_model), _bool2str(self.var_cmu), _bool2str(self.quadratic), _bool2str(self.cubic), _bool2str(self.scalar_nl), _bool2str(self.spec_nl), _bool2str(self.species), _bool2str(self.enthalpy)))        
        ncc_input_data.append(' '.join(['mean_on'.center(25), 'residual_smoothing_on'.center(25)]))
        ncc_input_data.append('{0:25s} {1:25s}'.format(_bool2str(self.mean_on), _bool2str(self.residual_smoothing_on)))  
        ncc_input_data.append(' '.join(['iter_turbulence'.center(25), 'iter_species'.center(25), 'iter_enthalpy'.center(25)]))
        ncc_input_data.append('{0:25s} {1:25s} {2:25s}'.format(_int2str(self.iter_turbulence), _int2str(self.iter_species), _int2str(self.iter_enthalpy)))   
        ncc_input_data.append(' '.join(['inv_sigma_k'.center(25), 'inv_sigma_e'.center(25)]))         
        ncc_input_data.append('{0:25.16f} {1:25.16f}'.format(self.inv_sigma_k, self.inv_sigma_e))  
        ncc_input_data.append(' '.join(['RCp'.center(25), 'C_mu'.center(25), 'C_e1'.center(25), 'C_e2'.center(25)]))        
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.RCp, self.C_mu, self.C_e1, self.C_e2))    
        ncc_input_data.append(' '.join(['laminar_viscosity'.center(25), 'effect_lam_visc'.center(25), 'under_relax'.center(25), 'RNG_cap'.center(25)]))        
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.laminar_viscosity, self.effect_lam_visc, self.under_relax, self.RNG_cap)) 
        ncc_input_data.append(' '.join(['1/Sch_1'.center(25), '1/Sch_2'.center(25), '1/Sch_3'.center(25), '1/Sch_4'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.inv_Sch_1, self.inv_Sch_2, self.inv_Sch_3, self.inv_Sch_4)) 
        ncc_input_data.append(' '.join(['stoic'.center(25), 'A1_EDC'.center(25), 'A2_EDC'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.stoic, self.A1_EDC, self.A2_EDC))         
        ncc_input_data.append(' '.join(['universal_gas_constant'.center(25), 'mw_fuel'.center(25), 'mw_oxidant'.center(25), 'mw_prod'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.universal_gas_constant, self.mw_fuel, self.mw_oxidant, self.mw_prod))         
        ncc_input_data.append(' '.join(['cp_fuel'.center(25), 'cp_oxidant'.center(25), 'cp_prod'.center(25), 'Prandtl#'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.cp_fuel, self.cp_oxidant, self.cp_prod, self.Prandtl))          
        ncc_input_data.append(' '.join(['combust'.center(25), 'chem_model'.center(25), 'href_fuel'.center(25), 'href_oxidant'.center(25), 'href_prod'.center(25), 'Tref'.center(25)]))
        ncc_input_data.append('{0:25s} {1:25s} {2:25.16f} {3:25.16f} {4:25.16f} {5:25.16f}'.format(_bool2str(self.combust), _enum2str(self.chem_model), self.href_fuel, self.href_oxidant, self.href_prod, self.Tref))         
        ncc_input_data.append(' '.join(['min_temperature'.center(25), 'max_temperature'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f}'.format(self.min_temperature, self.max_temperature))         
        ncc_input_data.append(' '.join(['lspray'.center(25), '# liquid injectors'.center(25)]))
        ncc_input_data.append('{0:25s} {1:25s}'.format(_bool2str(self.lspray), _int2str(self.num_liquids)))        
        ncc_input_data.append(' '.join(['2nd-order smoothing; aero'.center(25), 'k-e'.center(25), 'species'.center(25), 'enthalpy'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.aero2, self.k_e2, self.species2, self.enthalpy2))                  
        ncc_input_data.append(' '.join(['4th-order smoothing; aero'.center(25), 'k-e'.center(25), 'species'.center(25), 'enthalpy'.center(25)]))     
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f}'.format(self.aero4, self.k_e4, self.species4, self.enthalpy4))     
        ncc_input_data.append(' '.join(['4th-order switches; aero'.center(25), 'k-e'.center(25), 'species'.center(25), 'enthalpy'.center(25)]))    
        ncc_input_data.append('{0:25s} {1:25s} {2:25s} {3:25s}'.format(_bool2str(self.aero4_switch), _bool2str(self.k_e4_switch), _bool2str(self.species4_switch), _bool2str(self.enthalpy4_switch)))   
        ncc_input_data.append(' '.join(['k2'.center(25), 'beta_scale'.center(25)])) 
        ncc_input_data.append('{0:25.16f} {1:25.16f}'.format(self.k2, self.beta_scale))    
        ncc_input_data.append(' '.join(['expn'.center(25)]))
        ncc_input_data.append('{0:25.16f}'.format(self.expn)) 
        ncc_input_data.append(' '.join(['high_parallel'.center(25)]))
        ncc_input_data.append('{0:25s}'.format(_bool2str(self.high_parallel))) 
        ncc_input_data.append(' '.join(['chem_turb_model'.center(25), 'rxn_low_temp'.center(25)]))        
        ncc_input_data.append('{0:25s} {1:25.16f}'.format(_int2str(self.chem_turb_model), self.rxn_low_temp))    
        ncc_input_data.append(' '.join(['y_fuel_fuel'.center(25), 'y_o2_fuel'.center(25), 'y_co2_fuel'.center(25), 'y_h2o_fuel'.center(25), 'y_n2_fuel'.center(25)])) 
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f}'.format(self.y_fuel_fuel, self.y_o2_fuel, self.y_co2_fuel, self.y_h2o_fuel, self.y_n2_fuel)) 
        ncc_input_data.append(' '.join(['y_fuel_oxid'.center(25), 'y_o2_oxid'.center(25), 'y_co2_oxid'.center(25), 'y_h2o_oxid'.center(25), 'y_n2_oxid'.center(25)])) 
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f}'.format(self.y_fuel_oxid, self.y_o2_oxid, self.y_co2_oxid, self.y_h2o_oxid, self.y_n2_oxid)) 
        ncc_input_data.append(' '.join(['debug_res_on'.center(25), 'debug_res_start'.center(25), 'debug_res_print'.center(25)])) 
        ncc_input_data.append('{0:25s} {1:25s} {2:25s}'.format(_bool2str(self.debug_res_on), _int2str(self.debug_res_start), _int2str(self.debug_res_print)))  
        ncc_input_data.append(' '.join(['pg_on'.center(25)]))       
        ncc_input_data.append('{0:25s}'.format(_bool2str(self.pg_on))) 
        ncc_input_data.append(' '.join(['hp_eos'.center(25)]))           
        ncc_input_data.append('{0:25s}'.format(_int2str(self.hp_eos))) 
        ncc_input_data.append(' '.join(['VNN_solid'.center(25)]))        
        ncc_input_data.append('{0:25.16f}'.format(self.vnn_solid)) 
        ncc_input_data.append(' '.join(['angle_axi'.center(25)]))
        ncc_input_data.append('{0:25.16f}'.format(self.angle_axi))         
        ncc_input_data.append(' '.join(['int_chem_solver'.center(25), 'substeps'.center(25)]))        
        ncc_input_data.append('{0:25s} {1:25s}'.format(_int2str(self.int_chem_solver), _int2str(self.substeps))) 
        ncc_input_data.append(' '.join(['newton_temp'.center(25), 'tol_temp'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f}'.format(self.newton_temp, self.tol_temp)) 
        ncc_input_data.append(' '.join(['calc_lam_diff'.center(25), 'mass_diff_in_heat_flux'.center(25)]))    
        ncc_input_data.append('{0:25s} {1:25s}'.format(_bool2str(self.calc_lam_diff), _bool2str(self.mass_diff_in_heat_flux)))
        ncc_input_data.append(' '.join(['max_iter_tstatic'.center(25), 'newton_tstatic'.center(25), 'tol_tstatic'.center(25)]))
        ncc_input_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.max_iter_tstatic, self.newton_tstatic, self.tol_tstatic))
        ncc_input_data.append(' '.join(['over_write'.center(25), 'add_mut_switch'.center(25)]))    
        ncc_input_data.append('{0:25s} {1:25s}'.format(_bool2str(self.over_write), _bool2str(self.add_mut_switch)))            
        ncc_input_data.append(' '.join(['ce3'.center(25), 'add_ce3_term'.center(25)]))    
        ncc_input_data.append('{0:25.16f} {1:25s}'.format(self.over_write, _bool2str(self.add_ce3_term)))           
        ncc_input_data.append(' '.join(['enth_variance_on'.center(25), 'iter_enth_variance'.center(25)]))    
        ncc_input_data.append('{0:25s} {1:25s}'.format(_bool2str(self.enth_variance_on), _int2str(self.iter_enth_variance)))           
        ncc_input_data.append(' '.join(['tprime_intensity'.center(25), 'g_decrease_factor'.center(25)]))    
        ncc_input_data.append('{0:25.16f} {1:25.16f}'.format(self.tprime_intensity, self.g_decrease_factor))           
        ncc_input_data.append(' '.join(['apdf_tprime_small'.center(25), 'apdf_tprime_large'.center(25)]))    
        ncc_input_data.append('{0:25.16f} {1:25.16f}'.format(self.apdf_tprime_small, self.apdf_tprime_large))           
        ncc_input_data.append(' '.join(['apdf_simpson_eps'.center(25), 'apdf_simpson_jmax'.center(25)]))    
        ncc_input_data.append('{0:25.16f} {1:25s}'.format(self.apdf_simpson_eps, _int2str(self.apdf_simpson_jmax)))           
        ncc_input_data.append(' '.join(['thermal_wallfct'.center(25)]))    
        ncc_input_data.append('{0:25s}'.format(_bool2str(self.thermal_wallfct)))           
        ncc_input_data.append(' '.join(['WRITE_DB_CONTROL'.center(25)]))    
        ncc_input_data.append('{0:25s}'.format(_int2str(self.WRITE_DB_CONTROL)))
        
        # ---------------------------------------
        # --- Construct ncc_bcflags data list ---
        # ---------------------------------------   
        ncc_bcflags_data = ['1 \'inlet-mdot-norm\'']
        ncc_bcflags_data.append(' '.join(['mdot'.center(25), 'Tstatic'.center(25), 'turb'.center(25), 'Lmix'.center(25), 'Mach_no'.center(25)]))  
        ncc_bcflags_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f}'.format(self.bc1_mdot, self.bc1_Tstatic, self.bc1_turb, self.bc1_Lmix, self.bc1_Mach_no)) 
        ncc_bcflags_data.append(' '.join(['O2'.center(25), 'N2'.center(25)]))        
        ncc_bcflags_data.append('{0:25.14f} {1:25.14f}'.format(self.bc1_O2, self.bc1_N2))    
        ncc_bcflags_data.append('2 \'exit\'')
        ncc_bcflags_data.append(' '.join(['Pstatic'.center(25), 'Tstatic'.center(25), 'subsonic'.center(25)]))
        ncc_bcflags_data.append('{0:25.16f} {1:25.16f} {2:25s} '.format(self.bc2_Pstatic, self.bc2_Tstatic, _bool2str(self.bc2_subsonic)))        
        ncc_bcflags_data.append('3 \'periodic-a\'')     
        ncc_bcflags_data.append('4 \'periodic-b\'')
        # ncc_bcflags_data.append('5 \'symmetry\'')    
        ncc_bcflags_data.append('5 \'inject-wall-mdot-norm\'')
        ncc_bcflags_data.append(' '.join(['model'.center(25), 'Twall'.center(25), 'mdot'.center(25), 'Aflow'.center(25), 'Tstatic'.center(25), 'turb'.center(25), 'Lmix'.center(25)]))  
        ncc_bcflags_data.append('{0:25s} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f} {5:25.16f} {6:25.16f}'.format(_int2str(self.bc7_model), self.bc7_Twall, self.bc7_mdot, self.bc7_Aflow, self.bc7_Tstatic, self.bc7_turb, self.bc7_Lmix)) 
        ncc_bcflags_data.append(' '.join(['O2'.center(25), 'N2'.center(25)]))        
        ncc_bcflags_data.append('{0:25.14f} {1:25.14f}'.format(self.bc1_O2, self.bc1_N2))        

        # ----------------------------------------------
        # --- Construct ncc_material_flags data list ---
        # ----------------------------------------------
        ncc_material_flags_data = []
        for vol in range(1,21):
            ncc_material_flags_data.append(str(vol) + ' \'fluid\'')  
            ncc_material_flags_data.append(' '.join(['axisymmetric rotation'.center(25)]))
            ncc_material_flags_data.append('{0:25.16f}'.format(self.axisymmetric_rotation))         
            ncc_material_flags_data.append(' '.join(['initial u'.center(25), 'initial v'.center(25), 'initial w'.center(25), 'initial Pstatic'.center(25), 'initial k'.center(25), 'initial e'.center(25), 'initial Tstatic'.center(25)]))
            ncc_material_flags_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f} {5:25.16f} {6:25.16f}'.format(self.u_init, self.v_init, self.w_init, self.Pstatic_init, self.k_init, self.e_init, self.Tstatic_init)) 
            ncc_material_flags_data.append(' '.join(['O2'.center(25), 'N2'.center(25)]))        
            ncc_material_flags_data.append('{0:25.14f} {1:25.14f}'.format(self.O2_init, self.N2_init))

        # ------------------------------------------------
        # --- Construct ncc_liquid_solver.in data list ---
        # ------------------------------------------------
        ncc_liquid_solver_data = [' '.join(['ldread'.center(25), 'ispray_mod'.center(25), 'when_start_spray (1-no of liquids)'.center(25)])] 
        ncc_liquid_solver_data.append('{0:25s} {1:25s} {2:25s}'.format(_bool2str(self.ldread), _int2str(self.ispray_mod), (str(self.when_start_spray) + ' ') * self.num_injectors))
        ncc_liquid_solver_data.append(' '.join(['iswim'.center(25), 'wall_film_thickness'.center(25), 'offset'.center(25), 'coeff_of_restitution'.center(25)])) 
        ncc_liquid_solver_data.append('{0:25s} {1:25.16f} {2:25.16f} {3:25.16f}'.format(_int2str(self.iswim), self.wall_film_thickness, self.offset, self.coeff_of_restitution))        
        ncc_liquid_solver_data.append(' '.join(['dtml'.center(25), 'dtgl'.center(25), 'dtil'.center(25)])) 
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.dtml, self.dtgl, self.dtil))      
        ncc_liquid_solver_data.append(' '.join(['evap_calc_model'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_int2str(self.evap_calc_mod)))        
        ncc_liquid_solver_data.append(' '.join(['sklim_constant'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f}'.format(self.sklim_constant))    
        ncc_liquid_solver_data.append(' '.join(['use_constant_rholm'.center(25)]))        
        ncc_liquid_solver_data.append('{0:25s}'.format(_bool2str(self.use_constant_rholm)))    
        ncc_liquid_solver_data.append(' '.join(['tol_bmgs'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f}'.format(self.tol_bmgs))    
        ncc_liquid_solver_data.append(' '.join(['tol_xtold_lower'.center(25), 'tol_xtold_upper'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f}'.format(self.tol_xtold_lower, self.tol_xtold_upper))        
        ncc_liquid_solver_data.append(' '.join(['tol_ymk_lower'.center(25), 'tol_ymk_upper'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f}'.format(self.tol_ymk_lower, self.tol_ymk_upper))    
        ncc_liquid_solver_data.append(' '.join(['tol_xt_tdrop'.center(25), 'tol_xt_tboil'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f}'.format(self.tol_xt_tdrop, self.tol_xt_tboil))
        ncc_liquid_solver_data.append(' '.join(['force_evap_mass_factor'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f}'.format(self.force_evap_mass_factor))        
        ncc_liquid_solver_data.append(' '.join(['num_calls_no_evap'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_int2str(self.num_calls_no_evap)))        
        ncc_liquid_solver_data.append(' '.join(['evaporation_on'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_bool2str(self.evaporation_on)))            
        ncc_liquid_solver_data.append(' '.join(['i_interp_gas_part'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_int2str(self.i_interp_gas_part)))    
        ncc_liquid_solver_data.append(' '.join(['tol_prob_bprime_sprips'.center(25), 'chas_1'.center(25), 'chas_2'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.tol_prob_bprime_sprips, self.bprime_chas_1, self.bprime_chas_2))        
        ncc_liquid_solver_data.append(' '.join(['tol_prob_devu1_sprips'.center(25), 'chas_1'.center(25), 'chas_2'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.tol_prob_devu1_sprips, self.devu1_chas_1, self.devu1_chas_2))        
        ncc_liquid_solver_data.append(' '.join(['tol_prob_sklim_1'.center(25), 'sklim_2'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f}'.format(self.tol_prob_sklim_1, self.sklim_2))        
        ncc_liquid_solver_data.append(' '.join(['tol_prob_speed_chg_1'.center(25), 'chg_2'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f}'.format(self.tol_prob_speed_chg_1, self.chg_2))
        ncc_liquid_solver_data.append(' '.join(['relax_gas_src_mass'.center(25), 'moment'.center(25), 'energy'.center(25)]))
        ncc_liquid_solver_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.relax_gas_src_mass, self.moment, self.energy))        
        ncc_liquid_solver_data.append(' '.join(['i_relax_gas_src_spray'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_int2str(self.i_relax_gas_src_spray)))            
        ncc_liquid_solver_data.append(' '.join(['max_spray_int_count'.center(25), 'spray_int_interval'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s} {1:25s}'.format(_int2str(self.max_spray_int_count), _int2str(self.spray_int_interval)))        
        ncc_liquid_solver_data.append(' '.join(['i_update_icrit'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_int2str(self.i_update_icrit)))        
        ncc_liquid_solver_data.append(' '.join(['stochastic_spray'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s}'.format(_bool2str(self.stochastic_spray)))    
        ncc_liquid_solver_data.append(' '.join(['max_part_searches'.center(25), 'max_part_search_exchanges'.center(25)]))    
        ncc_liquid_solver_data.append('{0:25s} {1:25s}'.format(_int2str(self.max_part_searches), _int2str(self.max_part_search_exchanges)))    
        ncc_liquid_solver_data.append(' '.join(['limit_dtml'.center(25), 'ncells_per_dtml'.center(25)]))
        ncc_liquid_solver_data.append('{0:25s} {1:25.16f}'.format(_bool2str(self.limit_dtml), self.ncells_per_dtml))        
        ncc_liquid_solver_data.append(' '.join(['keyword_partsearch'.center(25)]))    
        ncc_liquid_solver_data.append(''.join(['{0:25s}'.format(('\''+ str(self.keyword_partsearch) + '\'').center(25))]))        
        ncc_liquid_solver_data.append(' '.join(['rklim'.center(25)]))    
        ncc_liquid_solver_data.append('{0:25.16f}'.format(self.rklim))

        # -------------------------------------------
        # --- Construct ncc_pm_input.in data list ---
        # -------------------------------------------
        ncc_pm_input_data = []
                
        if self.reacting or self.reacting2:
            ncc_pm_input_data = ['pbcfl = ' + str(self.pbcfl)]
            ncc_pm_input_data.append('etau_beta = ' + str(self.etau_beta))

        for location in ignition_loc:
             
            x = float(location[0])
            y = float(location[1])
            z = float(location[2])
                                              
            self.xmin1 = x + 0.0075
            self.xmax1 = self.xmin1 + 0.001
            self.ymin1 = y
            self.ymax1 = y
            self.zmin1 = z
            self.zmax1 = z

            if self.reacting:
                ncc_pm_input_data.append('glowplug=0,{0:.5f},{1:0.5f},{2:0.5f},{3:.5f},{4:0.5f},{5:.5f},{6:1s},{7:4s},{8:1s},{9:2s},{10:4s},-.007'.format(self.xmin1, self.ymin1, self.zmin1, self.xmax1, self.ymax1, self.zmax1, str(self.delT), str(self.T_off), str(self.ignite_off_model), str(self.ignite_on_restart), str(self.limit_ignition_duration)))

        ncc_pm_input_data.append('end')

        # ---------------------------------------------
        # --- Construct ncc_injector.in.# data list ---
        # ---------------------------------------------
        for index, injector in enumerate(injector_loc):
            
            self.x = float(injector[0]) + 0.0015 # --- Shift 1.5-mm downstream
            self.y = float(injector[1])
            self.z = float(injector[2])
  
            ncc_injector_data = [' '.join(['num_of_liquid_components'.center(25)])]  
            ncc_injector_data.append(' '.join(['Fuel Injector Input Info'.center(25)]))            
            ncc_injector_data.append('{0:25s}'.format(_int2str(self.number_of_liquid_components)))
            ncc_injector_data.append(' '.join(['symbols & mass fractions'.center(25)]))            
            ncc_injector_data.append('{0:25s}'.format(self.fuel_symbol.center(25)))
            ncc_injector_data.append('{0:25s}'.format(_int2str(1.00)))            
            ncc_injector_data.append(' '.join(['initial_liquid_temp'.center(25), 'superheat_vap_model'.center(25), 'Jet_A_properties'.center(25), 'init_mass_fraction_vapor'.center(25)]))    
            ncc_injector_data.append('{0:25.16f} {1:25s} {2:25s} {3:25.16f}'.format(self.initial_liquid_temperature, _int2str(self.superheat_vap_model), _bool2str(self.Jet_A_properties), self.init_mass_fraction_vapor))
            ncc_injector_data.append(' '.join(['atomization'.center(25), 'drop_breakup_model'.center(25), 'line_injection'.center(25), 'spray_table'.center(25), 'steady_spray_model'.center(25)]))
            ncc_injector_data.append('{0:25s} {1:25s} {2:25s} {3:25s} {4:25s}'.format(_bool2str(self.atomization), _bool2str(self.drop_breakup_model), _bool2str(self.line_injection), _bool2str(self.spray_table), _bool2str(self.steady_spray_model)))
            ncc_injector_data.append(' '.join(['no_of_holes'.center(25), 'no_of_streams'.center(25), 'no_of_droplet_groups'.center(25)]))
            ncc_injector_data.append('{0:25s} {1:25s} {2:25s}'.format(_int2str(self.no_of_holes), _int2str(self.no_of_streams), _int2str(self.no_of_droplet_groups)))
            ncc_injector_data.append(' '.join(['lmdis'.center(25), 'smdm'.center(25), '3D_cone'.center(25), 'size_min'.center(25), 'size_max'.center(25), 'stochastic_injection'.center(25)]))            
            ncc_injector_data.append('{0:25s} {1:25.16f} {2:25s} {3:25.16f} {4:25.16f} {5:25s}'.format(_int2str(self.lmdis), self.smdm, _bool2str(self.cone), self.size_min, self.size_max, _bool2str(self.stochastic_injection)))
            ncc_injector_data.append(' '.join(['x'.center(25), 'y'.center(25), 'z'.center(25), 'mdot_spray'.center(25), 'v_inj'.center(25), 'alpha'.center(25), 'beta'.center(25), 'cone angle'.center(25), 'delta cone angle'.center(25), 'swrl_angle'.center(25)]))    
            ncc_injector_data.append('{0:25.16f} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f} {5:25.16f} {6:25.16f} {7:25.16f} {8:25.16f} {9:25.16f}'.format(self.x, self.y, self.z, self.mdot_spray, self.v_inj, self.alpha, self.beta, self.theta, self.dtheta, self.swrl_angle)) 
            ncc_injector_data.append(' '.join(['atom_type'.center(25), 'breakup_type'.center(25), 'dia_hole'.center(25), 'delp_inj'.center(25), 'liq_vel'.center(25), 'gas_u'.center(25), 'gas_v'.center(25), 'gas_w'.center(25), 'pcl_start'.center(25), 'pcl_end'.center(25)]))    
            ncc_injector_data.append('{0:25s} {1:25s} {2:25.16f} {3:25.16f} {4:25.16f} {5:25.16f} {6:25.16f} {7:25.16f} {8:25s} {9:25s}'.format(_int2str(self.atom_type), _int2str(self.breakup_type), self.dia_hole, self.delp_inj, self.liq_vel, self.gas_u, self.gas_v, self.gas_w, _int2str(self.pcl_start), _int2str(self.pcl_end)))                        
            ncc_injector_data.append(' '.join(['line_injection_type'.center(25), 'xi'.center(25), 'yi'.center(25), 'zi'.center(25), 'xf'.center(25), 'yf'.center(25), 'zf'.center(25), 'ri'.center(25)]))           
            ncc_injector_data.append('{0:25s} {1:25.16f} {2:25.16f} {3:25.16f} {4:25.16f} {5:25.16f} {6:25.16f} {7:25.16f}'.format(_int2str(self.line_injection_type), self.xi, self.yi, self.zi, self.xf, self.yf, self.zf, self.ri))               
            
            
            # -----------------------------------------        
            # --- Write output to ncc_injector.in.x ---
            # -----------------------------------------
            filename = Path('OpenMDAO') + '/ncc_injector.in.' + str(index + 1)
            file_handle = open(filename, 'w')
            file_handle.write("\n".join(ncc_injector_data))
            file_handle.close()        
 
        # -----------------------------------------
        # --- Construct nodal results data list ---
        # -----------------------------------------   
        nodal_results_data = ['number of files'.center(25)]
        nodal_results_data.append(' '.join(['1'.center(25)]))
        nodal_results_data.append(' '.join(['--------------------------'.center(25)]))
        nodal_results_data.append(' '.join(['filename 1'.center(25)]))
        nodal_results_data.append('{0:25s}'.format(self.results_fn + '.dat'))
        nodal_results_data.append(' '.join(['number of variables'.center(25)]))
        nodal_results_data.append(' '.join(['16'.center(25)]))
        nodal_results_data.append(' '.join(['variable names'.center(25)]))
        nodal_results_data.append("'u' 'v' 'w' 'Vmag' 'M' 'pg' 'p' 'Ptot_approx' 'rho' 'T' 'Ttot_approx' 'k' 'e' 'mu' 'Unmixedness' 'yall'")             

        # ----------------------------------------------
        # --- Construct ncc_results_aux.in data list ---
        # ----------------------------------------------   
        ncc_results_aux_data = ['write_aux'.center(25)]
        ncc_results_aux_data.append(' '.join(['T'.center(25)])) 
        ncc_results_aux_data.append(' '.join(['number_aux_vars'.center(25)]))
        ncc_results_aux_data.append(' '.join(['1'.center(25)]))
        ncc_results_aux_data.append(' '.join(['parameter names'.center(25)]))
        
        if self.reacting2:
            ncc_results_aux_data.append('{0:25s}'.format("'unmixed' 'C11H21'".center(25))) 
        else:
            ncc_results_aux_data.append('{0:25s}'.format("'unmixed' 'C12H23'".center(25)))            

        # ------------------------------------------------
        # --- Construct ncc_walltime_stop.in data list ---
        # ------------------------------------------------      
        ncc_walltime_stop_data = [' '.join(['twall_run (hrs)'.center(25), 'twall_restart_files (mins)'.center(25), 'twall_factor (<=1)'.center(25)])]
        ncc_walltime_stop_data.append('{0:25.16f} {1:25.16f} {2:25.16f}'.format(self.twall_run, self.twall_restart_files, self.twall_factor))  
     
        # -------------------------------------------------
        # --- Construct ncc_boundary_write.in data list ---
        # -------------------------------------------------
        ncc_boundary_data = ['Write boundary face data (T or F)'.center(50)]
        ncc_boundary_data.append(' '.join(['T'.center(25)]))        
        ncc_boundary_data.append(' '.join(['Num BC flags below'.center(25)]))
        ncc_boundary_data.append(' '.join(['5'.center(25)]))        
        ncc_boundary_data.append(' '.join(['BC-flag'.center(25), 'no_of_parameters'.center(25), 'parameters_to_print'.center(40)]))   
        ncc_boundary_data.append(' '.join(['1'.center(25), '2'.center(25), "'htot_avg' 'ptot_avg'".center(40)]))
        ncc_boundary_data.append(' '.join(['3'.center(25), '1'.center(25), "'htot_avg'".center(40)]))
        ncc_boundary_data.append(' '.join(['4'.center(25), '1'.center(25), "'htot_avg'".center(40)]))        
        ncc_boundary_data.append(' '.join(['5'.center(25), '2'.center(25), "'ptot_avg' 'T_avg'".center(40)]))

        if self.reacting2:
            ncc_boundary_data.append(' '.join(['2'.center(25), '5'.center(25), "'htot_avg' 'ptot_avg' 'T_avg' 'EI_NOx' 'EI_CO'".center(40)]))
            ncc_boundary_data.append('{0:25.16f}'.format(self.mdot_spray * self.num_injectors))            
        else:
            ncc_boundary_data.append(' '.join(['2'.center(25), '3'.center(25), "'htot_avg' 'ptot_avg' 'T_avg'".center(40)]))
             
        # --------------------------------------
        # --- Construct post-1d.in data list ---
        # --------------------------------------
        post1d_data = ['Geometry'.center(25)]
        post1d_data.append(' '.join(['3d'.center(25)]))        
        post1d_data.append(' '.join(['Num BC flags below'.center(25)]))
        post1d_data.append(' '.join(['4'.center(25)]))        
        post1d_data.append(' '.join(['BC-flag'.center(25), 'no_of_parameters'.center(25), 'parameters_to_print'.center(40)]))   
        post1d_data.append(' '.join(['1'.center(25), '2'.center(25), "'htot_avg' 'ptot_avg'".center(40)]))
        post1d_data.append(' '.join(['3'.center(25), '1'.center(25), "'htot_avg'".center(40)]))
        post1d_data.append(' '.join(['4'.center(25), '1'.center(25), "'htot_avg'".center(40)]))        
        
        if self.reacting2:
            post1d_data.append(' '.join(['2'.center(25), '5'.center(25), "'htot_avg' 'ptot_avg' 'T_avg' 'EI_NOx' 'EI_CO'".center(40)]))              
        else:
            post1d_data.append(' '.join(['2'.center(25), '3'.center(25), "'htot_avg' 'ptot_avg' 'T_avg'".center(40)]))
            
        post1d_data.append('{0:25.16f}'.format(self.mdot_spray * self.num_injectors))            
     
        # ==================================================================================================================
        # ===                                   Write NCC input data lists to file                                       ===
        # ==================================================================================================================     
        # ------------------------------------
        # --- Write output to ncc_input.in ---
        # ------------------------------------          
        filename = Path('OpenMDAO') + '/ncc_input.in'    
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_input_data))
        file_handle.close()
        
        # ----------------------------------------        
        # --- Write output to ncc_chemistry.in ---
        # ----------------------------------------
        filename = Path('OpenMDAO') + '/ncc_chemistry.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_chemistry_data))
        file_handle.close()
        
        # --------------------------------------        
        # --- Write output to ncc_bcflags.in ---
        # --------------------------------------
        filename = Path('OpenMDAO') + '/ncc_bcflags.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_bcflags_data))
        file_handle.close()
        
        # ---------------------------------------------        
        # --- Write output to ncc_material_flags.in ---
        # ---------------------------------------------
        filename = Path('OpenMDAO') + '/ncc_material_flags.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_material_flags_data))
        file_handle.close()
        
        # --------------------------------------------       
        # --- Write output to ncc_liquid_solver.in ---
        # --------------------------------------------
        filename = Path('OpenMDAO') + '/ncc_liquid_solver.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_liquid_solver_data))
        file_handle.close()         
        
        # --------------------------------------------     
        # --- Write output to ncc_ignition.in ---
        # --------------------------------------------
        # filename = Path('OpenMDAO') + '/ncc_ignition.in'
        # file_handle = open(filename, 'w')
        # file_handle.write("\n".join(ncc_ignition_data))
        # file_handle.close() 

        # --------------------------------------------     
        # --- Write output to ncc_pm_input.in ---
        # --------------------------------------------
        filename = Path('OpenMDAO') + '/ncc_pm_input.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_pm_input_data))
        file_handle.close() 
            
        # ----------------------------------------        
        # --- Write output to nodal_results.in ---
        # ----------------------------------------
        filename = Path('OpenMDAO') + '/nodal_results.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(nodal_results_data))
        file_handle.close()         
        
        # ------------------------------------------        
        # --- Write output to ncc_results_aux.in ---
        # ------------------------------------------
        filename = Path('OpenMDAO') + '/ncc_results_aux.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_results_aux_data))
        file_handle.close() 
        
        # ------------------------------------------        
        # --- Write output to ncc_walltime_stop.in ---
        # ------------------------------------------
        filename = Path('OpenMDAO') + '/ncc_walltime_stop.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_walltime_stop_data))
        file_handle.close()   

        # ---------------------------------------------        
        # --- Write output to ncc_boundary_write.in ---
        # ---------------------------------------------
        filename = Path('OpenMDAO') + '/ncc_boundary_write.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(ncc_boundary_data))
        file_handle.close()
        
        # ------------------------------------------        
        # --- Write output to post-1d.in ---
        # ------------------------------------------
        filename = Path('OpenMDAO') + '/post-1d.in'
        file_handle = open(filename, 'w')
        file_handle.write("\n".join(post1d_data))
        file_handle.close()        

        if self._clean_start:
            shutil.copy2(Path('Patran') + '\Config' + self.config.split('_')[0] + '/patran.out', Path('OpenMDAO'))
            
        # ------------------      
        # --- Submit Job ---
        # ------------------
        print 'Starting Job on NAS ...'             

        self.external_files = [FileMetadata('ncc_input.in', input=True),                                   
                               FileMetadata('ncc_bcflags.in', input=True),
                               FileMetadata('ncc_id_period.in', input=True),                            
                               FileMetadata('ncc_chemistry.in', input=True),
                               FileMetadata('ncc_material_flags.in', input=True),                                   
                               FileMetadata('nodal_results.in', input=True),
                               FileMetadata('ncc_boundary_write.in', input=True),
                               FileMetadata('open_ncc_boundary.out.*', output=True),                                
                               FileMetadata('post-1d.in', input=True),
                               FileMetadata('post-1d.out', output=True),                                
                               FileMetadata('ncc_walltime_stop.in', input=True), 
                               # FileMetadata('ncc_results_aux.in', input=True),                                                        Uncomment this line when the aux vars are reintroduced
                               FileMetadata('qsub.py', input=True),
                               FileMetadata('preprocess', input=True),                                       
                               FileMetadata('postprocess', input=True),
                               FileMetadata('write_macro.py', input=True),                                   
                               FileMetadata('tecplot.dat', output=True),                   
                               FileMetadata('open_ncc_resd.out', output=True),
                               FileMetadata('*.png', output=True, binary=True),                                  
                               FileMetadata('Ts.txt', output=True),                                   
                               FileMetadata('TKE.txt', output=True)]                  
        
        if self.reacting:            
            self.external_files.extend([FileMetadata('ncc_injector.in.*', input=True),
                                        FileMetadata('ncc_liquid_solver.in', input=True), 
                                        FileMetadata('ncc_pm_input.in', input=True),          
                                        FileMetadata('run_reacting', input=True),                                                       
                                        FileMetadata('tecplot_unsteady_droplets.dat', output=True)])
        
        elif self.reacting2:
            self.external_files.extend([FileMetadata('ncc_injector.in.*', input=True),
                                        FileMetadata('ncc_liquid_solver.in', input=True),
                                        FileMetadata('ncc_pm_input.in', input=True),
                                        FileMetadata('ncc_converge.in', input=True),                                         
                                        FileMetadata('ncc_converge.out', output=True),              
                                        FileMetadata('run_reacting2', input=True),              
                                        FileMetadata('tecplot_unsteady_droplets.dat', output=True),
                                        FileMetadata('change_results_species.in', input=True),
                                        FileMetadata('Unmixedness.txt', output=True),
                                        FileMetadata('NO.txt', output=True)])                  
        
        else:
            self.external_files.extend([FileMetadata('run_nonreacting', input=True)])                                   
            
        if self._clean_start:
            shutil.copy2(Path('Patran') + '\Config' + self.config.split('_')[0] + '/patran.out', Path('OpenMDAO'))
            self.external_files.extend([FileMetadata('patran.out', input=True)])              

        self.resources = {'localhost': False}
            
        # -----------------------------
        # --- Execute the component ---
        # -----------------------------
        super(NCCcomponent, self).execute()  

        print 'Job Complete ...'
        
        self._clean_start = False        
        
if __name__ == "__main__":
 
    enable_console() 
    logging.getLogger().setLevel(logging.DEBUG) 
  
    # --- Create 4 allocators (reducancy for reliability) adding each to the manager ---
    allocator1 = NAS_Allocator(name='PFE20-DMZ1',
                              dmz_host='dmzfs1.nas.nasa.gov',
                              server_host='pfe20')
                                          
    allocator2 = NAS_Allocator(name='PFE20-DMZ2',
                              dmz_host='dmzfs2.nas.nasa.gov',
                              server_host='pfe20')   

    allocator3 = NAS_Allocator(name='PFE21-DMZ1',
                              dmz_host='dmzfs1.nas.nasa.gov',
                              server_host='pfe21')  

    allocator4 = NAS_Allocator(name='PFE21-DMZ2',
                              dmz_host='dmzfs2.nas.nasa.gov',
                              server_host='pfe21')                              
                              
    RAM.add_allocator(allocator1)
    RAM.add_allocator(allocator2)
    RAM.add_allocator(allocator3)
    RAM.add_allocator(allocator4)    

    NCC_Comp = NCCcomponent()
    NCC_Comp.max_iterations_per_time_step = 100
    NCC_Comp.nonreacting = True
    
    try:
        NCC_Comp.run()
    finally:
        # Without this a process is left polling at the server.
        allocator.shutdown()   
