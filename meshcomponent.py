''' OpenMDAO mesh component '''
''' Updates geometry via the SolidWorks API and exports as a STEP file''' 
''' Loads STEP file into cubit and applies a geometry adaptive tetrahedral mesh'''

# --- Inherent python/system level imports
import os

# --- External python library imports (i.e. matplotlib, numpy, scipy)
import cubit
import shutil
import math

# --- OpenMDAO imports
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Str, Int, Float
from pathname import Path

cubit.init(" ")

def run_cubit(config, num_injectors, cut_plane):
    ''' Function containing python commands to run Cubit and generate a ~1,000,000 element tetrahedral mesh '''
    
    # --- Create New Cubit File
    cubit.cmd('new')

    # --- Import File with Heal Option specified
    fname = Path('STEP') + '\Config' + config + '.step'
    cubit.cmd('import step "%s" heal'%fname)

    cubit.cmd('Set Max Memory On 200')  
    cubit.cmd('regularize volume all')
    cubit.cmd('undo on')    
    cubit.cmd('set multisweep off')
    if num_injectors == 1:
        cubit.cmd('webcut volume all with plane yplane offset ' + str(cut_plane - 0.5) + ' imprint merge ')
        cubit.cmd('webcut volume 2 with plane yplane offset 6.5 imprint merge ')
        cubit.cmd('volume all scheme tetmesh')
        cubit.cmd('volume 2 sizing function type skeleton scale 3.5 time_accuracy_level 2.5 min_size auto max_size 0.08 max_gradient 1.125 min_depth 3 max_depth 6 min_num_layers_3d 3 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 2')
        cubit.cmd('volume 1 sizing function type skeleton scale 3.0 time_accuracy_level 2.5 min_size 0.01 max_size 0.225 max_gradient 1.1 min_depth 2.0 max_depth 7 min_num_layers_3d 1 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 1')        
        cubit.cmd('volume 3 sizing function type skeleton scale 3.0 time_accuracy_level 2.5 min_size 0.01 max_size 0.225 max_gradient 1.1 min_depth 2.0 max_depth 7 min_num_layers_3d 1 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 3')         
    elif num_injectors == 7:
        cubit.cmd('webcut volume all with plane yplane offset ' + str(cut_plane - 0.4) + ' imprint merge ')
        cubit.cmd('webcut volume 2 with plane yplane offset 6.0 imprint merge ')
        cubit.cmd('volume all scheme tetmesh')
        cubit.cmd('volume 2 sizing function type skeleton scale 5.0 time_accuracy_level 2.5 min_size auto max_size 0.06 max_gradient 1.2 min_depth 3 max_depth 6 min_num_layers_3d 4 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 2')  
        cubit.cmd('volume 1 sizing function type skeleton scale 3.0 time_accuracy_level 2.5 min_size 0.02 max_size 0.2 max_gradient 1.125 min_depth 2.0 max_depth 7 min_num_layers_3d 1 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 1')        
        cubit.cmd('volume 3 sizing function type skeleton scale 3.0 time_accuracy_level 2.5 min_size 0.02 max_size 0.2 max_gradient 1.125 min_depth 2.0 max_depth 7 min_num_layers_3d 1 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 3')        
    else:
        cubit.cmd('webcut volume all with plane yplane offset ' + str(cut_plane - 0.3) + ' imprint merge ')
        cubit.cmd('webcut volume 2 with plane yplane offset 5.5 imprint merge ')    
        cubit.cmd('volume all scheme tetmesh')
        cubit.cmd('volume 2 sizing function type skeleton scale 5.0 time_accuracy_level 2.5 min_size auto max_size 0.04 max_gradient 1.2 min_depth 3 max_depth 6 min_num_layers_3d 4 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 2')  
        cubit.cmd('volume 1 sizing function type skeleton scale 3.0 time_accuracy_level 2.5 min_size 0.01 max_size 0.175 max_gradient 1.1 min_depth 2.0 max_depth 7 min_num_layers_3d 1 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 1')        
        cubit.cmd('volume 3 sizing function type skeleton scale 3.0 time_accuracy_level 2.5 min_size 0.01 max_size 0.175 max_gradient 1.1 min_depth 2.0 max_depth 7 min_num_layers_3d 1 min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')
        cubit.cmd('mesh volume 3') 

    cubit.cmd('undo group begin')
    cubit.cmd('volume all smooth scheme condition number beta 1.5 cpu 1.0')    
    cubit.cmd('smooth volume all')
    cubit.cmd('undo group end')
    cubit.cmd('create pressure 1 on surface ( at 8.881784e-016 0.000000e+000 0.000000e+000 ordinal 1 ordered ) magnitude 1')
    cubit.cmd('create pressure 2 on surface ( at 2.775558e-017 1.40000e+001 -1.072052e-034 ordinal 1 ordered ) magnitude 2')

    cubit.cmd('save as "' + Path('Cubit') + '\Config' + config + '.cub" overwrite')
    cubit.cmd('export patran "' + Path('OpenMDAO') + '\patran.out" overwrite')
    cubit.cmd('reset')

    

class MeshComp(Component):
    ''' OpenMDAO component for meshing geometry '''  

    config = Str('30', iotype ='in', desc ='current design configuration number')    
    num_injectors = Int(7, iotype = 'in', desc = 'number of fuel injector modules')
    injector_dia = Float(0.8, iotype = 'in', desc = 'injector diameter')
    venturi_angle = Float(40, iotype = 'in', desc = 'angle of converging diverging venturi')    
    
    def __init__(self, *args, **kwargs):
        # ------------------------------------------
        # --- Constructor for the mesh component ---
        # ------------------------------------------
        super(MeshComp, self).__init__(*args, **kwargs)        
        self.force_execute = True
           
    def execute(self):
    
        print 'Updating Local Directory Structure ...'        
        shutil.rmtree(Path('Patran') + '\Config' + self.config, True)
        os.mkdir(Path('Patran') + '\Config' + self.config)
        
        print 'Generating Mesh for Configuration ' + self.config + ' ...'
        
        dome_location = 14.0 - 9.5
        half_venturi_length = (self.injector_dia * 0.5 - self.injector_dia * 0.32)/math.tan(self.venturi_angle * math.pi / 180.0)
        throat_thickness = self.injector_dia * 0.04125 
        vane_length = self.injector_dia * 0.52125
        cut_plane = dome_location - 2 * half_venturi_length - throat_thickness - vane_length
        
        # --- Executes Cubit commands to generate an unstructured tetrahedral mesh ---
        run_cubit(self.config, self.num_injectors, cut_plane)  

        shutil.copy2(Path('OpenMDAO') + '/patran.out', Path('Patran') + '\Config' + self.config[0])        

if __name__ == "__main__":
    
    #--- Default test case ---      
    Mesh_Comp = MeshComp()
    Mesh_Comp.run()

   