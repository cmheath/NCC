from openmdao.main.api import create_io_traits, set_as_top
from openmdao.lib.datatypes.api import Bool
from npsscomponent import NPSScomponent
from pathname import Path
from math import pi

class NpssReacting(NPSScomponent):

    def __init__(self):

        super(NpssReacting, self).__init__(arglist='NPSS_Combustor_Model.run')       
        
    def configure(self):
        super(NpssReacting, self).configure()
        self.run_command = 'run_reacting()'
        self.force_execute = True

        create_io_traits(self, ['Inlet.Pt',
                                'Inlet.Tt',
                                'Inlet.W',
                                'Inlet.Fl_O.Aphy',
                                'Conv_Duct.Fl_O.Aphy',
                                'Div_Duct.Fl_O.Aphy',
                                'Fuel.Tfuel',
                                'Comb.Fl_O.Aphy',                                
                                'Comb.dPqPBase',
                                'Comb.effBase',
                                'Bld1.BldOut.fracW',
                                'Bld1.Fl_O.Aphy', 
                                'Bld2.Fl_O.Aphy',
                                'Globals.T4',
                                'Globals.FAR'
                               ], iotype='in')

        create_io_traits(self, ['Inlet.Fl_O.V', 
                                'Inlet.Fl_O.Ts',
                                'Inlet.Fl_O.W',
                                'Inlet.Fl_O.rhos',
                                'Conv_Duct.Fl_O.W',                                
                                'Comb.Fl_O.Ps',
                                'Comb.Fl_O.Ts',                                
                                'Comb.FAR',
                                'Bld1.BldOut.W',
                                'Bld1.Fl_O.Ts'                                
                               ], iotype='out')
                               
class NpssNonreacting(NPSScomponent):

    def __init__(self):

        super(NpssNonreacting, self).__init__(arglist='NPSS_Combustor_Model.run')

    def configure(self):
        super(NpssNonreacting, self).configure()
        self.run_command = 'run_nonreacting()'
        self.force_execute = True
        
        create_io_traits(self, ['Inlet.Pt',
                                'Inlet.Tt',
                                'Inlet.W',
                                'Inlet.Fl_O.Aphy',
                                'Conv_Duct.Fl_O.Aphy',
                                'Div_Duct.Fl_O.Aphy',
                                'Fuel.Tfuel',
                                'Comb.Fl_O.Aphy',                                
                                'Comb.dPqPBase',
                                'Comb.effBase',
                                'Bld1.BldOut.fracW',
                                'Bld1.Fl_O.Aphy', 
                                'Bld2.Fl_O.Aphy',
                                'Globals.T4'                                
                               ], iotype='in')

        create_io_traits(self, ['Inlet.Fl_O.V', 
                                'Inlet.Fl_O.Ts',
                                'Inlet.Fl_O.W',
                                'Inlet.Fl_O.rhos',
                                'Conv_Duct.Fl_O.W',
                                'Comb.Fl_O.Ps',
                                'Comb.Fl_O.Ts',                                
                                'Comb.FAR',
                                'Bld1.BldOut.W',
                                'Bld1.Fl_O.Ts'                                
                               ], iotype='out')                          

if __name__ == '__main__':

    reacting = set_as_top(NpssReacting())
    
    venturi_throat_dia = 0.55                    
    venturi_throat_area = pi*(venturi_throat_dia/2)**2 * 3
    half_sector_angle = 7.5
    inner_dia = 16.5
    outer_dia = 21.5
    sector_area = (pi*(outer_dia/2)**2 - pi*(inner_dia/2)**2) * half_sector_angle*2/360.0
    injector_area = pi*(1.0/2)**2 * 2 + pi*(1.3/2)**2
    
    m_dot = 45.99 / 24.0 # 24 sectors in a full annular combustor
    
    reacting.Inlet_Fl_O_Aphy = sector_area
    reacting.Conv_Duct_Fl_O_Aphy = venturi_throat_area
    reacting.Div_Duct_Fl_O_Aphy = injector_area
    reacting.Comb_Fl_O_Aphy = sector_area
    reacting.Bld1_Fl_O_Aphy = sector_area    
    reacting.Bld2_Fl_O_Aphy = sector_area 
    
    reacting.Globals_T4 = 3200
    reacting.Inlet_W = m_dot
    reacting.Bld1_BldOut_fracW = 0.2
    reacting.run()

