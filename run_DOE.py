''' Top level assembly for the Lean Direct Injection (LDI) design space analysis '''

# --- Inherent python/system level imports
import os
import sys
import logging
import numpy as np
from math import pi

# --- OpenMDAO main and library imports
from openmdao.main.api import Assembly, SequentialWorkflow, enable_console
from openmdao.main.resource import ResourceAllocationManager as RAM
from openmdao.lib.drivers.api import DOEdriver, FixedPointIterator, IterateUntil, CaseIteratorDriver
from openmdao.lib.casehandlers.api import CSVCaseIterator, DumpCaseRecorder
from openmdao.lib.doegenerators.api import OptLatinHypercube
from nas_access import NAS_Allocator

# --- OpenMDAO component imports
from countercomponent import DOECounterComp, SubCounterComp
from restartcomponent import RestartComp
from geometrycomponent import GeometryComp
from remotemeshcomponent import MeshComp
from npsscomp import NpssReacting, NpssNonreacting
from ncccomponent import NCCcomponent
from tecplotcomponent import TecplotComp
from pathname import Path
from init_npss import init_dP
from run_local import run_local

class Analysis(Assembly):
    ''' Top level for the Lean Direct Injection (LDI) combustion design space analysis '''

    def configure(self): 

        # --------------------------------------------------------------------------- #
        # --- Instantiate Counter Component
        # --------------------------------------------------------------------------- #          
        self.add('counter', DOECounterComp())
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate Subcounter Component
        # --------------------------------------------------------------------------- #          
        self.add('subcounter', SubCounterComp()) 
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate Restart Component
        # --------------------------------------------------------------------------- #          
        self.add('restart_doe', RestartComp()) 
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate NPSS Non-Reacting Component
        # --------------------------------------------------------------------------- #       
        self.add('npss_nonreacting', NpssNonreacting())
        self.npss_nonreacting.Comb_dPqPBase = init_dP()

        # --------------------------------------------------------------------------- #        
        # --- Instantiate NPSS Reacting Component
        # --------------------------------------------------------------------------- #  
        self.add('npss_reacting', NpssReacting())       
        self.npss_reacting.Comb_dPqPBase = init_dP()

        # --------------------------------------------------------------------------- #
        # --- Instantiate Geometry Component
        # --------------------------------------------------------------------------- #  
        self.add('geometry', GeometryComp()) 

        # --------------------------------------------------------------------------- #
        # --- Instantiate Mesh Component
        # --------------------------------------------------------------------------- #        
        self.add('mesh', MeshComp())

        # --------------------------------------------------------------------------- #
        # --- Instantiate NCC Non-reacting Component (1-Step Chemistry)
        # --------------------------------------------------------------------------- #  
        self.add('ncc_nonreacting', NCCcomponent())
        self.ncc_nonreacting.max_iterations_per_time_step = 40000
        self.ncc_nonreacting.nonreacting = True  
        self.ncc_nonreacting.continuity_goal = 2000    
        self.ncc_nonreacting.restarts = 40000        
        self.ncc_nonreacting.CFL = 0.8
        self.ncc_nonreacting.mass_imbalance_goal = 1.0E-2
        self.ncc_nonreacting.aero2 = -1.0E-3
        self.ncc_nonreacting.k_e2 = -1.0E-3
        self.ncc_nonreacting.species2 = -1.0E-3
        self.ncc_nonreacting.enthalpy2 = -1.0E-3
        self.ncc_nonreacting.aero4 = 0.05
        self.ncc_nonreacting.k_e4 = 0.05
        self.ncc_nonreacting.species4 = 0.05
        self.ncc_nonreacting.enthalpy4 = 0.05
        self.ncc_nonreacting.twall_run = 4.0
        self.ncc_nonreacting._num_procs = 400 # --- Must be divisible by 20 to run on Ivy Bridge
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate NCC Reacting Component (1-Step Chemistry)
        # --------------------------------------------------------------------------- #  
        self.add('ncc_reacting', NCCcomponent())
        self.ncc_reacting.reacting = True
        self.ncc_reacting.continuity_goal = 750
        self.ncc_reacting.max_iterations_per_time_step = 15000
        self.ncc_reacting.restarts = 15000
        self.ncc_reacting.CFL = 0.8
        self.ncc_reacting.pbcfl = 0.375
        self.ncc_reacting.etau_beta = 0.15        
        self.ncc_reacting.mass_imbalance_goal = 1.0E-3
        self.ncc_reacting.aux_var_name_list = "'unmixed' 'C12H23'"
        self.ncc_reacting.spec_name_ignite = 'C12H23'
        self.ncc_reacting.fuel_symbol = 'C12H23'
        self.ncc_reacting.combust = True
        self.ncc_reacting.lspray = True
        self.ncc_reacting.aero2 = -2.5E-3
        self.ncc_reacting.k_e2 = -2.5E-3
        self.ncc_reacting.species2 = -2.5E-3
        self.ncc_reacting.enthalpy2 = -2.5E-3
        self.ncc_reacting.aero4 = 0.05
        self.ncc_reacting.k_e4 = 0.05
        self.ncc_reacting.species4 = 0.05
        self.ncc_reacting.enthalpy4 = 0.05
        self.ncc_reacting.ignite_done_model = 0
        self.ncc_reacting.ignition_on = True
        self.ncc_reacting.no_of_streams = 16
        self.ncc_reacting.twall_run = 4.0
        self.ncc_reacting._num_procs = 400 # --- Must be divisible by 20 to run on Ivy Bridge
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate NCC Reacting Component (Detailed Chemistry)
        # --------------------------------------------------------------------------- #  
        self.add('ncc_reacting2', NCCcomponent())        
        self.ncc_reacting2.reacting2 = True
        self.ncc_reacting2.continuity_goal = 750
        self.ncc_reacting2.max_iterations_per_time_step = 10000
        self.ncc_reacting2.restarts = 10000
        self.ncc_reacting2.CFL = 0.9
        self.ncc_reacting.pbcfl = 0.5
        self.ncc_reacting2.mass_imbalance_goal = 1.0E-3
        self.ncc_reacting2.aux_var_name_list = "'unmixed' 'C11H21'" 
        self.ncc_reacting2.spec_name_ignite = 'C11H21'
        self.ncc_reacting2.fuel_symbol = 'C11H21'
        self.ncc_reacting2.combust = True
        self.ncc_reacting2.lspray = True
        self.ncc_reacting2.aero2 = -1.0E-3
        self.ncc_reacting2.k_e2 = -1.0E-3
        self.ncc_reacting2.species2 = -1.0E-3
        self.ncc_reacting2.enthalpy2 = -1.0E-3
        self.ncc_reacting2.aero4 = 0.05
        self.ncc_reacting2.k_e4 = 0.05
        self.ncc_reacting2.species4 = 0.05
        self.ncc_reacting2.enthalpy4 = 0.05
        self.ncc_reacting2.energy = 0.0001
        self.ncc_reacting2.when_start_spray = 1
        self.ncc_reacting2.ignite_done_model = 0
        self.ncc_reacting2.ignition_on = True
        self.ncc_reacting2.twall_run = 7.5
        self.ncc_reacting2._num_procs = 400 # --- Must be divisible by 20 to Run on Ivy Bridge
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate Tecplot Nonreacting Component
        # --------------------------------------------------------------------------- #  
        self.add('tecplot_nonreacting', TecplotComp())
        self.tecplot_nonreacting.nonreacting = True
        
        # --------------------------------------------------------------------------- #
        # --- Instantiate Tecplot Reacting Spray Component
        # --------------------------------------------------------------------------- #  
        self.add('tecplot_reacting', TecplotComp())        
        self.tecplot_reacting.reacting = True

        # --------------------------------------------------------------------------- #
        # --- Instantiate Tecplot Reacting Spray Component
        # --------------------------------------------------------------------------- #  
        self.add('tecplot_reacting2', TecplotComp())        
        self.tecplot_reacting2.reacting2 = True
        
        # --------------------------------------------------------------------------- #
        # --- Create Driver Instances
        # --------------------------------------------------------------------------- #       
        # --- Top Level Assembly Driver
        self.add('driver', IterateUntil())
        self.driver.max_iterations = 1
        self.driver.workflow = SequentialWorkflow()
       
        # --- Inner Nonreacting Driver (Fixed Point)      
        self.add('nonreacting_driver', FixedPointIterator())
        self.nonreacting_driver.workflow = SequentialWorkflow()
        self.nonreacting_driver.step_size = 0.125
        self.nonreacting_driver.max_iteration = 1
        self.nonreacting_driver.tolerance = 0.001

        # --- Inner Reacting Driver
        self.add('reacting_driver', IterateUntil())
        self.reacting_driver.max_iterations = 1
        self.reacting_driver.workflow = SequentialWorkflow()
        
        # --- Inner Reacting Driver #2
        self.add('reacting_driver2', IterateUntil())
        self.reacting_driver2.max_iterations = 2
        self.reacting_driver2.workflow = SequentialWorkflow()
        self.reacting_driver2.add_stop_condition('tecplot_reacting2.dTqTBase < 0.0025')
        
        # --- Run Design at All ICAO Power Settings
        self.add('power_hook', CaseIteratorDriver())
        self.power_hook.workflow = SequentialWorkflow()        
        self.power_hook.iterator = CSVCaseIterator(filename='ICAOsettings.csv')
        self.power_hook.workflow.add(['reacting_driver2'])
        
        # --------------------------------------------------------------------------- #
        # --- Create Main Assembly Workflow
        # --------------------------------------------------------------------------- #        
        # --- Add component instances to top-level assembly
        #self.driver.workflow.add(['counter', 'restart_doe', 'geometry', 'mesh', 'nonreacting_driver', 'reacting_driver', 'reacting_driver2'])
        #self.driver.workflow.add(['counter', 'restart_doe', 'geometry', 'nonreacting_driver', 'reacting_driver', 'reacting_driver2'])
        #self.driver.workflow.add(['counter', 'restart_doe', 'geometry', 'reacting_driver2'])
        #self.driver.workflow.add(['counter', 'restart_doe', 'geometry', 'npss_nonreacting', 'npss_reacting'])
        self.driver.workflow.add(['counter', 'restart_doe', 'geometry', 'power_hook'])
        
        # --------------------------------------------------------------------------- #
        # --- Create Sub-Assembly Workflows
        # --------------------------------------------------------------------------- #
        # --- Inner Nonreacting Loop - Solve via fixed point iteration
        self.nonreacting_driver.workflow.add(['subcounter', 'npss_nonreacting', 'ncc_nonreacting', 'tecplot_nonreacting'])
    
        # --- Add solver independents and dependents for fixed point iterator
        self.nonreacting_driver.add_parameter(['npss_nonreacting.Comb_dPqPBase', 'npss_reacting.Comb_dPqPBase'], low = -9.e99, high = 9.e99)      
        self.nonreacting_driver.add_constraint('tecplot_nonreacting.dPqPBase = npss_nonreacting.Comb_dPqPBase')
        
        # --- Inner Reacting Loop - Run Once
        self.reacting_driver.workflow.add(['subcounter', 'npss_reacting', 'ncc_reacting', 'tecplot_reacting'])

        # --- Inner Reacting Loop #2 - Run Once
        self.reacting_driver2.workflow.add(['subcounter', 'npss_reacting', 'ncc_reacting2', 'tecplot_reacting2'])
   
        # --------------------------------------------------------------------------- #
        # --- Add Driver Events --- Comment this section out to override numbering
        # --------------------------------------------------------------------------- #
        if (self.subcounter._iteration == 0):
            self.driver.add_event('subcounter.reset_iteration')
            self.driver.add_event('ncc_nonreacting.clean_start')
            
            self.power_hook.add_event('subcounter.reset_iteration')
            self.power_hook.add_event('ncc_nonreacting.clean_start')
        
        # --------------------------------------------------------------------------- #        
        # --- Specify Case Recorders
        # --------------------------------------------------------------------------- #            
        self.driver.case_outputs = ['geometry.pilot_recession', 'geometry.vane_height', 'geometry.venturi_angle', 'nonreacting_driver.ncc_nonreacting.FAR', 'nonreacting_driver.tecplot_nonreacting.dPqPBase']  
        
        self.power_hook.recorders = [DumpCaseRecorder()]
                
        # --------------------------------------------------------------------------- #        
        # --- Create Data Connections 
        # --------------------------------------------------------------------------- #
        self.connect('counter.config', ['subcounter.case', 'restart_doe.config', 'geometry.config', 'mesh.config'])
        self.connect('subcounter.config', ['ncc_nonreacting.config', 'ncc_reacting.config', 'ncc_reacting2.config', 'tecplot_nonreacting.config', 'tecplot_reacting.config', 'tecplot_reacting2.config'])        
        self.connect('restart_doe.pilot_recession', ['geometry.pilot_recession', 'mesh.pilot_recession'])        
        self.connect('restart_doe.venturi_angle', ['geometry.venturi_angle', 'mesh.venturi_angle'])
        self.connect('restart_doe.vane_height', ['geometry.vane_height', 'mesh.vane_height'])        
        
        self.connect('geometry.injector_dia',  ['ncc_nonreacting.injector_dia', 'ncc_nonreacting.bc1_Lmix', 'ncc_reacting.injector_dia', 'ncc_reacting.bc1_Lmix', 'ncc_reacting2.injector_dia', 'ncc_reacting2.bc1_Lmix'])
        self.connect('geometry.sector_area', ['npss_nonreacting.Inlet_Fl_O_Aphy', 'npss_nonreacting.Comb_Fl_O_Aphy', 'npss_nonreacting.Bld1_Fl_O_Aphy', 'npss_nonreacting.Bld2_Fl_O_Aphy'])
        self.connect('geometry.sector_area', ['npss_reacting.Inlet_Fl_O_Aphy', 'npss_reacting.Comb_Fl_O_Aphy', 'npss_reacting.Bld1_Fl_O_Aphy', 'npss_reacting.Bld2_Fl_O_Aphy'])        
        self.connect('geometry.venturi_throat_area', ['npss_nonreacting.Conv_Duct_Fl_O_Aphy', 'npss_reacting.Conv_Duct_Fl_O_Aphy'])
        self.connect('geometry.venturi_exit_area', ['npss_nonreacting.Div_Duct_Fl_O_Aphy', 'npss_reacting.Div_Duct_Fl_O_Aphy'])      
        
        self.connect('npss_nonreacting.Conv_Duct_Fl_O_W', 'ncc_nonreacting.bc1_mdot') 
        self.connect('npss_nonreacting.Inlet_Fl_O_V', 'ncc_nonreacting.u_init')   
        self.connect('npss_nonreacting.Inlet_Fl_O_Ts', ['ncc_nonreacting.bc1_Tstatic', 'ncc_nonreacting.Tstatic_init'])       
        self.connect('npss_nonreacting.Comb_Fl_O_Ps', ['ncc_nonreacting.bc2_Pstatic', 'ncc_nonreacting.Pstatic_init']) 
        self.connect('npss_nonreacting.Comb_Fl_O_Ts', 'ncc_nonreacting.bc2_Tstatic') 
        self.connect('npss_nonreacting.Comb_FAR', 'ncc_nonreacting.FAR')          
        self.connect('npss_nonreacting.Inlet_Fl_O_rhos', 'ncc_nonreacting.rho_air')
        self.connect('npss_nonreacting.Bld1_BldOut_W', 'ncc_nonreacting.bc7_mdot')
        self.connect('npss_nonreacting.Bld1_Fl_O_Ts', 'ncc_nonreacting.bc7_Tstatic')        
        
        self.connect('npss_reacting.Conv_Duct_Fl_O_W', ['ncc_reacting.bc1_mdot', 'ncc_reacting2.bc1_mdot'])
        self.connect('npss_reacting.Inlet_Fl_O_Ts', ['ncc_reacting.bc1_Tstatic', 'ncc_reacting2.bc1_Tstatic'])       
        self.connect('npss_reacting.Comb_Fl_O_Ps', ['ncc_reacting.bc2_Pstatic', 'ncc_reacting2.bc2_Pstatic', 'ncc_reacting.Pstatic_init', 'ncc_reacting2.Pstatic_init']) 
        self.connect('npss_reacting.Comb_Fl_O_Ts', ['ncc_reacting.bc2_Tstatic', 'ncc_reacting2.bc2_Tstatic', 'ncc_reacting.Tstatic_init', 'ncc_reacting2.Tstatic_init'])
        self.connect('npss_reacting.Comb_FAR', ['ncc_reacting.FAR', 'ncc_reacting2.FAR'])                          
        self.connect('npss_reacting.Inlet_Fl_O_rhos', ['ncc_reacting.rho_air', 'ncc_reacting2.rho_air'])
        self.connect('npss_reacting.Bld1_BldOut_W', ['ncc_reacting.bc7_mdot', 'ncc_reacting2.bc7_mdot'])        
        self.connect('npss_reacting.Bld1_Fl_O_Ts', ['ncc_reacting.bc7_Tstatic', 'ncc_reacting2.bc7_Tstatic'])
        
if __name__ == '__main__':
  
    enable_console()
    logging.getLogger().setLevel(logging.DEBUG)

    # --------------------------------------------------------------------------- #        
    # --- Execute OpenMDAO DOE Component ---
    # --------------------------------------------------------------------------- #          
    top_level_analysis = Analysis()    
    
    try:
        top_level_analysis.run()
    finally:
        print 'allocator shutdown'   

    # --------------------------------END---------------------------------------- #        
       