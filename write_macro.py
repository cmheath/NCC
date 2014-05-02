# --- Inherent python/system level imports
import os
from sys import argv

# --- External python library imports (i.e. matplotlib, numpy, scipy)
import numpy as np

def _TKEcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to TKE. '''
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("kt") + 1))  # --- FIX-ME - Wrong Variable
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 1800, 10)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')     
        
    return tecplot_macro

def _Vcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to V. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("V<sub>Axial</sub> (m/s)") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(-10, 80, 10)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro
    
def _Tcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to Tstatic. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("T<sub>Static</sub> (K)") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(900, 2500, 17)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro    
    
def _Unmixedcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to Unmixedness. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("mu") + 1))  # --- FIX-ME - Wrong Variable
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 1, 11)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    

    return tecplot_macro  

def _Fuelcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to Fuel. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("Jet-A") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 0.5, 11)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro   
    
def _NOcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to NO. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("NO") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0.000, 0.00012, 13)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro 

def _COcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to CO. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("CO") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 0.12, 13)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro

def _CO2contour(tecplot_macro, varlist):
    ''' Sets plot contour variable to CO2. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("CO2") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 0.25, 11)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro

def _OHcontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to NO. '''    
    tecplot_macro.append('$!GLOBALCONTOUR 1  VAR = %d \n' % (varlist.index("OH") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 1 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 0.0025, 11)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    
    
    return tecplot_macro    
    
def _Diametercontour(tecplot_macro, varlist):
    ''' Sets plot contour variable to Diameter. '''
    tecplot_macro.append('$!GLOBALCONTOUR 2  VAR = %d \n' % (varlist.index("DIA") + 1))
    tecplot_macro.append('$!CONTOURLEVELS NEW \n') 
    tecplot_macro.append('CONTOURGROUP = 2 \n') 
    tecplot_macro.append('RAWDATA \n') 
    contour_levels = np.linspace(0, 20.0, 11)
    
    tecplot_macro.append(str(len(contour_levels)) + ' \n')
    
    for contour_val in contour_levels:
        tecplot_macro.append(str(contour_val) + ' \n')    

    return tecplot_macro  

    ''' OpenMDAO component for Processing Flow Solution Data in Tecplot '''
    
def WriteMacro(state):

    nonreacting = False
    reacting = False
    reacting2 = False
    
    if (state == '/run_nonreacting'):
        nonreacting = True
        print 'nonreacting'
    elif (state == '/run_reacting'):
        reacting = True
        print 'reacting'
    else:
        reacting2 = True
        print 'reacting2'
    
    # Decompose Slice Regions to Cylindrical Coordinates
    root = os.getcwd()

    axial = [-0.01, -0.005, 1E-4, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035] # axial locations at which planar slices are extracted
 
    if reacting or nonreacting:
        #varlist = ["X", "Y", "Z", "u", "v", "w", "Vmag", "M", "pg", "p", "Ptot_approx", "rho", "T", "Ttot_approx", "k", "e", "mu", "Unmixedness", "C12H23", "O2", "CO2", "H2O", "N2"]
        varlist = ["X", "Y", "Z", "V<sub>Axial</sub> (m/s)", "v", "w", "Pg", "turb_k", "turb_e", "u_tau", "T<sub>Static</sub> (K)", "rho", "enthalpy", "Mach-No.", "mu", "cp", "kt", "Jet-A", "O2", "CO2", "H2O", "N2"]
        varlist_init = ["X", "Y", "Z", "u", "v", "w", "Pg", "turb_k", "turb_e", "u_tau", "t", "rho", "enthalpy", "Mach-No.", "mu", "cp", "kt", "C12H23", "O2", "CO2", "H2O", "N2"]        
    else:
        varlist = ["X", "Y", "Z", "V<sub>Axial</sub> (m/s)", "v", "w", "Pg", "turb_k", "turb_e", "u_tau", "T<sub>Static</sub> (K)", "rho", "enthalpy", "Mach-No.", "mu", "cp", "kt", "Jet-A", "CO2", "CO", "O2", "O", "OH", "N2", "H2", "H", "H2O", "HO2", "NO", "N2O", "N", "CH"]
        varlist_init = ["X", "Y", "Z", "u", "v", "w", "Pg", "turb_k", "turb_e", "u_tau", "t", "rho", "enthalpy", "Mach-No.", "mu", "cp", "kt", "C11H21", "CO2", "CO", "O2", "O", "OH", "N2", "H2", "H", "H2O", "HO2", "NO", "N2O", "N", "CH"]
        
    varlist2 = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"]
    
    # ------------------------------------------------
    # --- Initialize values for ProcessResults.mcr ---
    # ------------------------------------------------  
    tecplot_macro = ['#!MC 1300\n']
    tecplot_macro.append('# Created by Tecplot 360 build 13.1.0.15185 \n')   
    tecplot_macro.append('$!VarSet |MFBD| = \'\\grcfsb\mydocs$\cheath1\' \n')

    tecplot_macro.append('$!READDATASET \'"tecplot.dat"\' \n')
    tecplot_macro.append('READDATAOPTION = NEW \n')
    tecplot_macro.append('RESETSTYLE = YES \n')
    tecplot_macro.append('INCLUDETEXT = NO \n')        
    tecplot_macro.append('INCLUDEGEOM = YES \n')
    tecplot_macro.append('INCLUDECUSTOMLABELS = NO \n')
    tecplot_macro.append('VARLOADMODE = BYNAME \n')
    tecplot_macro.append('ASSIGNSTRANDIDS = YES \n')
    tecplot_macro.append('INITIALPLOTTYPE = CARTESIAN3D \n')
    tecplot_macro.append('VARNAMELIST = \'"' + '", "'.join(varlist_init) + '"\' \n')
    tecplot_macro.append('$!FIELDMAP [1]  SURFACES{SURFACESTOPLOT = EXPOSEDCELLFACES} \n')        
    tecplot_macro.append('$!FIELDLAYERS SHOWSHADE = YES \n')
    tecplot_macro.append('$!FIELDLAYERS USELIGHTINGEFFECT = YES \n')
    tecplot_macro.append('$!FIELDLAYERS USETRANSLUCENCY = YES \n')
    tecplot_macro.append('$!FIELDMAP [1]  EFFECTS{SURFACETRANSLUCENCY = 75} \n')
    tecplot_macro.append('$!FIELDMAP [1]  EFFECTS{LIGHTINGEFFECT = GOURAUD} \n')
    tecplot_macro.append('$!GLOBALEDGE MINCREASEANGLE = 155 \n')
    
    tecplot_macro.append('$!RENAMEDATASETVAR \n')
    tecplot_macro.append('  VAR = 4 \n') 
    tecplot_macro.append('  NAME = \'V<sub>Axial</sub> (m/s)\' \n')
    tecplot_macro.append('$!RENAMEDATASETVAR \n')    
    tecplot_macro.append('  VAR = 11 \n') 
    tecplot_macro.append('  NAME = \'T<sub>Static</sub> (K)\' \n') 
    tecplot_macro.append('$!RENAMEDATASETVAR \n')    
    tecplot_macro.append('  VAR = 18 \n') 
    tecplot_macro.append('  NAME = \'Jet-A\' \n')
    
    if not nonreacting and os.path.isfile('tecplot_unsteady_droplets.dat'):
        tecplot_macro.append('$!READDATASET \'"tecplot_unsteady_droplets.dat"\' \n')
        tecplot_macro.append('READDATAOPTION = APPEND \n')
        tecplot_macro.append('RESETSTYLE = NO \n')
        tecplot_macro.append('INCLUDETEXT = NO \n')        
        tecplot_macro.append('INCLUDEGEOM = NO \n')
        tecplot_macro.append('INCLUDECUSTOMLABELS = NO \n')
        tecplot_macro.append('VARLOADMODE = BYNAME \n')
        tecplot_macro.append('ASSIGNSTRANDIDS = YES \n')
        tecplot_macro.append('INITIALPLOTTYPE = CARTESIAN3D \n')            
        varlist.append("DIA")
        varlist.append("N_drop")
        varlist.append("Liquid")
        tecplot_macro.append('VARNAMELIST = \'"' + '", "'.join(varlist) + '"\' \n')                        
        tecplot_macro.append('$!FIELDMAP [1]  SCATTER{SHOW = NO} \n')
        tecplot_macro.append('$!GLOBALSCATTER VAR = 4 \n')
        tecplot_macro.append('$!FIELDMAP [2]  SCATTER{SYMBOLSHAPE{GEOMSHAPE = SPHERE}} \n')
        tecplot_macro.append('$!FIELDLAYERS SHOWSCATTER = YES \n')
        
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{X = 0.057150002568960856} \n')
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Y = 0.24259006703844532} \n')
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Z = 0.14577820474321776} \n')
        tecplot_macro.append('$!THREEDVIEW PSIANGLE = 0 \n')
        tecplot_macro.append('$!THREEDVIEW THETAANGLE = 0 \n')
        tecplot_macro.append('$!THREEDVIEW ALPHAANGLE = 0 \n')        
        tecplot_macro.append('$!VIEW DATAFIT \n')         
        tecplot_macro.append('$!VIEW SETMAGNIFICATION \n')
        tecplot_macro.append('MAGNIFICATION = 0.8 \n')
        tecplot_macro.append('$!SET3DEYEDISTANCE \n')
        tecplot_macro.append('  EYEDISTANCE = 15 \n')        

        tecplot_macro.append('$!THREEDAXIS XDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS XDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')         
        tecplot_macro = _Diametercontour(tecplot_macro, varlist)
        tecplot_macro.append('$!FIELDMAP [2]  SCATTER{FILLCOLOR = MULTI2} \n')
        tecplot_macro.append('$!FIELDMAP [2]  SCATTER{FRAMESIZE = 0.3} \n')  
        tecplot_macro.append('$!FIELDMAP [2]  SCATTER{SIZEBYVARIABLE = NO} \n')
        tecplot_macro.append('$!FIELDMAP [2]  SCATTER{COLOR = MULTI2} \n')            
        tecplot_macro.append('$!REDRAW \n')            
        tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
        tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'SideSprayDistribution.png') + '\' \n')
        tecplot_macro.append('$!EXPORT \n')
        tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')

        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{X = 0.85274513616587977} \n')        
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Y = -0.0059070000000000919} \n')
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Z = 0.24043000000000012} \n')
        tecplot_macro.append('$!THREEDVIEW PSIANGLE = 90 \n')        
        tecplot_macro.append('$!THREEDVIEW THETAANGLE = -90 \n')            
        tecplot_macro.append('$!THREEDVIEW ALPHAANGLE = 0 \n')
        tecplot_macro.append('$!ROTATE3DVIEW X \n')
        tecplot_macro.append('  ANGLE = 90 \n')
        tecplot_macro.append('  ROTATEORIGINLOCATION = DEFINEDORIGIN \n')

        tecplot_macro.append('$!VIEW DATAFIT \n')            
        tecplot_macro.append('$!VIEW SETMAGNIFICATION \n')
        tecplot_macro.append('MAGNIFICATION = 0.7 \n')
        tecplot_macro.append('$!REDRAW \n')     
        tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
        tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'EndSprayDistribution.png') + '\' \n')
        tecplot_macro.append('$!EXPORT \n')
        tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n') 
        tecplot_macro.append('$!FIELDLAYERS SHOWSCATTER = NO \n')        
        tecplot_macro.append('$!DELETEZONES [2] \n')      
    
    for theta in [1, 2]:
    
        # --- Turn on TKE Contour          
        tecplot_macro = _TKEcontour(tecplot_macro, varlist)
        tecplot_macro.append('$!FIELDLAYERS SHOWCONTOUR = YES \n')
        tecplot_macro.append('$!FIELDLAYERS USELIGHTINGEFFECT = NO \n')
        tecplot_macro.append('$!FIELDLAYERS USETRANSLUCENCY = NO \n')        
        tecplot_macro.append('$!FIELDLAYERS SHOWEDGE = NO \n')
        
        if (theta == 1):
            # Create Pilot Side-View Slice from Plane and Orient View
            tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{X = 0}} \n')
            tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{Z = 1}} \n')
            tecplot_macro.append('$!CREATESLICEZONEFROMPLANE \n')
            tecplot_macro.append('SLICESOURCE = VOLUMEZONES \n')
            tecplot_macro.append('FORCEEXTRACTIONTOSINGLEZONE = YES \n')
            tecplot_macro.append('COPYCELLCENTEREDVALUES = NO \n')
            tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{X = 0.057150002568960856} \n')
            tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Y = 0.24259006703844532} \n')
            tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Z = 0.14577820474321776} \n')
            tecplot_macro.append('$!THREEDVIEW PSIANGLE = 0 \n')
            tecplot_macro.append('$!THREEDVIEW THETAANGLE = 0 \n')
            tecplot_macro.append('$!THREEDVIEW ALPHAANGLE = 0 \n')
            fname = 'PilotView'
            
        else:        
            # Create Main Side-View Slice from Plane    
            tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{Y = -0.10279253678724681}} \n')
            tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{Z = 0.99470281711717423}} \n')
            
            tecplot_macro.append('$!CREATESLICEZONEFROMPLANE \n')
            tecplot_macro.append('SLICESOURCE = VOLUMEZONES \n')                
            tecplot_macro.append('FORCEEXTRACTIONTOSINGLEZONE = YES \n')
            tecplot_macro.append('COPYCELLCENTEREDVALUES = NO \n')
            tecplot_macro.append('$!ROTATE3DVIEW X \n')
            tecplot_macro.append('  ANGLE = -6.0 \n')
            tecplot_macro.append('  ROTATEORIGINLOCATION = DEFINEDORIGIN \n')
            fname = 'MainView'       
                        
        tecplot_macro.append('$!ACTIVEFIELDMAPS = [2] \n')                   
        tecplot_macro.append('$!VIEW DATAFIT \n')
        tecplot_macro.append('$!VIEW SETMAGNIFICATION \n')
        tecplot_macro.append('MAGNIFICATION = 0.9 \n')        
        tecplot_macro.append('$!SET3DEYEDISTANCE \n')
        tecplot_macro.append('  EYEDISTANCE = 5 \n')
        tecplot_macro.append('$!THREEDAXIS XDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS XDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')        
        tecplot_macro.append('$!THREEDAXIS XDETAIL{AUTOGRID = NO} \n')
        tecplot_macro.append('$!THREEDAXIS XDETAIL{GRSPACING = 0.01} \n')
        tecplot_macro.append('$!VIEW CENTER \n')        
        tecplot_macro.append('$!REDRAW \n')     
        tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
        tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'TKE.png') + '\' \n')
        tecplot_macro.append('$!EXPORT \n')
        tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')  

        # --- Turn on V Contour and Save Image
        tecplot_macro = _Vcontour(tecplot_macro, varlist)            
        tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
        tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'Vel.png') + '\' \n')
        tecplot_macro.append('$!EXPORT \n')
        tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
        
        if reacting or reacting2:
            # --- Turn on T Contour and Save Image
            tecplot_macro = _Tcontour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'Tmp.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')        
     
            # --- Turn on Unmixedness Contour and Save Image 
            tecplot_macro = _Unmixedcontour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'Unmix.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')

            # --- Turn on Fuel Contour and Save Image
            tecplot_macro = _Fuelcontour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'Fuel.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')    

        if reacting2:
            # --- Turn on NO Contour and Save Image
            tecplot_macro = _NOcontour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'NO.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')   

            # --- Turn on CO Contour and Save Image
            tecplot_macro = _COcontour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'CO.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')

            # --- Turn on CO2 Contour and Save Image
            tecplot_macro = _CO2contour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'CO2.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
            
            # --- Turn on OH Contour and Save Image
            tecplot_macro = _OHcontour(tecplot_macro, varlist)            
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')        
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1199 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, fname + 'OH.png') + '\' \n')
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')            

        tecplot_macro.append('$!DELETEZONES [2] \n')
        tecplot_macro.append('$!ACTIVEFIELDMAPS = [1] \n')        

    tecplot_macro = _TKEcontour(tecplot_macro, varlist)
    
    tecplot_macro.append('$!FIELDLAYERS SHOWEDGE = YES \n')    
    tecplot_macro.append('$!FIELDLAYERS USELIGHTINGEFFECT = NO \n')
    tecplot_macro.append('$!FIELDLAYERS USETRANSLUCENCY = YES \n')
    tecplot_macro.append('$!FIELDMAP [1]  CONTOUR{SHOW = NO} \n')
    
    tecplot_macro.append('$!THREEDAXIS YDETAIL{SHOWAXIS = YES} \n')
    tecplot_macro.append('$!THREEDAXIS ZDETAIL{SHOWAXIS = YES} \n')
    tecplot_macro.append('$!WORKSPACEVIEW FITPAPER \n')        

    for X1 in axial:        
        #--- Export Axial Slice Image Data
        tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{X = 500}} \n')
        tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{Y = 0}} \n')
        tecplot_macro.append('$!GLOBALTHREED SLICE{ORIGIN{X = ' + str(X1) + '}} \n')
        tecplot_macro.append('$!CREATESLICEZONEFROMPLANE  \n')
        tecplot_macro.append('SLICESOURCE = VOLUMEZONES \n')
        tecplot_macro.append('FORCEEXTRACTIONTOSINGLEZONE = YES \n')
        tecplot_macro.append('COPYCELLCENTEREDVALUES = NO \n')
        
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{X = 0.85274513616587977} \n')        
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Y = -0.0059070000000000919} \n')
        tecplot_macro.append('$!THREEDVIEW VIEWERPOSITION{Z = 0.24043000000000012} \n')
        tecplot_macro.append('$!THREEDVIEW PSIANGLE = 90 \n')        
        tecplot_macro.append('$!THREEDVIEW THETAANGLE = -90 \n')            
        tecplot_macro.append('$!THREEDVIEW ALPHAANGLE = 0 \n')
        tecplot_macro.append('$!ROTATE3DVIEW X \n')
        tecplot_macro.append('  ANGLE = 90 \n')
        tecplot_macro.append('  ROTATEORIGINLOCATION = DEFINEDORIGIN \n') 
        tecplot_macro.append('$!THREEDAXIS YDETAIL{GRIDLINES{SHOW = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{GRIDLINES{SHOW = NO}} \n')  

        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TICKS{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TICKLABEL{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS ZDETAIL{TITLE{SHOWONAXISLINE = NO}} \n')
        tecplot_macro.append('$!THREEDAXIS YDETAIL{TITLE{SHOWONAXISLINE = NO}} \n')

        tecplot_macro.append('$!FIELDMAP [2]  EFFECTS{SURFACETRANSLUCENCY = 12} \n')
        tecplot_macro.append('$!FIELDMAP [1]  EFFECTS{SURFACETRANSLUCENCY = 99} \n')
        tecplot_macro.append('$!GLOBALCONTOUR 1  LEGEND{SHOW = YES} \n')
        tecplot_macro.append('$!VIEW DATAFIT \n')
        tecplot_macro.append('$!VIEW SETMAGNIFICATION \n')
        tecplot_macro.append('MAGNIFICATION = 0.75 \n')
        tecplot_macro.append('$!VIEW TRANSLATE \n')
        tecplot_macro.append('X = -10 \n')
        tecplot_macro.append('Y = -10 \n')
        tecplot_macro.append('$!SET3DEYEDISTANCE \n')
        tecplot_macro.append('  EYEDISTANCE = 0.1 \n')        
        tecplot_macro.append('$!REDRAWALL  \n')            
        tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
        tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'TKE_%02d.png') % round(X1 * 1000) + '\' \n') 
        tecplot_macro.append('$!EXPORT \n')
        tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')            
        tecplot_macro = _Vcontour(tecplot_macro, varlist)
        tecplot_macro.append('$!REDRAWALL  \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
        tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
        tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'Vel_%02d.png') % round(X1 * 1000) + '\' \n')
        tecplot_macro.append('$!EXPORT \n')
        tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
        
        if reacting or reacting2:
            tecplot_macro = _Tcontour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'Tmp_%02d.png') % round(X1 * 1000) + '\' \n') 
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')            
            tecplot_macro = _Unmixedcontour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n') 
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'Unmix_%02d.png') % round(X1 * 1000) + '\' \n') 
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
            tecplot_macro = _Fuelcontour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n') 
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'Fuel_%02d.png') % round(X1 * 1000) + '\' \n')  
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')  
        
        if reacting2:
            tecplot_macro = _NOcontour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n') 
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'NO_%02d.png') % round(X1 * 1000) + '\' \n')  
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
            tecplot_macro = _COcontour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n') 
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'CO_%02d.png') % round(X1 * 1000) + '\' \n')  
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
            tecplot_macro = _CO2contour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n') 
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'CO2_%02d.png') % round(X1 * 1000) + '\' \n')  
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')
            tecplot_macro = _OHcontour(tecplot_macro, varlist)
            tecplot_macro.append('$!REDRAWALL  \n') 
            tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
            tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 1769 \n')
            tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'OH_%02d.png') % round(X1 * 1000) + '\' \n')  
            tecplot_macro.append('$!EXPORT \n')
            tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n')                
        
        tecplot_macro = _TKEcontour(tecplot_macro, varlist)
        tecplot_macro.append('$!REDRAWALL  \n')                        
        tecplot_macro.append('$!DELETEZONES [2] \n')   

    for X in np.linspace(0.0, 0.19, 191):
        tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{X = 500}} \n')
        tecplot_macro.append('$!GLOBALTHREED SLICE{NORMAL{Y = 0}} \n')        
        tecplot_macro.append('$!GLOBALTHREED SLICE{ORIGIN{X = ' + str(X) + '}} \n')
        tecplot_macro.append('$!CREATESLICEZONEFROMPLANE \n')
        tecplot_macro.append('SLICESOURCE = VOLUMEZONES \n')
        tecplot_macro.append('FORCEEXTRACTIONTOSINGLEZONE = YES \n')
        tecplot_macro.append('COPYCELLCENTEREDVALUES = NO \n')

    # FIX-ME - Variable Indices Wrong
    tecplot_macro.append('$!EXTENDEDCOMMAND \n') 
    tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")
    tecplot_macro.append(r"COMMAND = 'Integrate [2-192] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=11 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'" + ' \n')        
    tecplot_macro.append('$!EXTENDEDCOMMAND  \n')
    tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")
    tecplot_macro.append("COMMAND = 'SaveIntegrationResults FileName = \"Ts.txt\"' \n")
    tecplot_macro.append('$!EXTENDEDCOMMAND  \n')
    tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")
    tecplot_macro.append(r"COMMAND = 'Integrate [2-192] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=15 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'" + ' \n')
    tecplot_macro.append('$!EXTENDEDCOMMAND \n')
    tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")
    tecplot_macro.append("COMMAND = 'SaveIntegrationResults FileName = \"TKE.txt\"' \n")
    tecplot_macro.append('$!EXTENDEDCOMMAND  \n')
    tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")
    tecplot_macro.append(r"COMMAND = 'Integrate [2-192] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=18 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'" + ' \n')
    tecplot_macro.append('$!EXTENDEDCOMMAND \n')
    tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")
    tecplot_macro.append("COMMAND = 'SaveIntegrationResults FileName = \"Unmixedness.txt\"' \n")      
    
    if reacting2:
        tecplot_macro.append('$!EXTENDEDCOMMAND  \n')
        tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")           
        tecplot_macro.append(r"COMMAND = 'Integrate [2-192] VariableOption=\'Average\' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=26 Absolute=\'F\' ExcludeBlanked=\'F\' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver=\'Cells\' IntegrateBy=\'Zones\' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults=\'F\' PlotAs=\'Result\' TimeMin=0 TimeMax=0'" + ' \n')
        tecplot_macro.append('$!EXTENDEDCOMMAND \n')
        tecplot_macro.append("COMMANDPROCESSORID = 'CFDAnalyzer3' \n")           
        tecplot_macro.append("COMMAND = 'SaveIntegrationResults FileName = \"NO.txt\"' \n")        
    
    tecplot_macro.append('$!DELETEZONES [2-192] \n')
      
    # --- Load and draw convergence history plot
    tecplot_macro.append('$!NEWLAYOUT \n')     
    tecplot_macro.append('$!READDATASET \'"open_ncc_resd.out"\' \n')
    tecplot_macro.append('READDATAOPTION = NEW \n')
    tecplot_macro.append('RESETSTYLE = YES \n')
    tecplot_macro.append('INCLUDETEXT = NO \n')
    tecplot_macro.append('INCLUDEGEOM = NO \n')
    tecplot_macro.append('INCLUDECUSTOMLABELS = NO \n')
    tecplot_macro.append('VARLOADMODE = BYNAME \n')
    tecplot_macro.append('ASSIGNSTRANDIDS = YES \n')
    tecplot_macro.append('INITIALPLOTTYPE = XYLINE \n')
    tecplot_macro.append('VARNAMELIST = \'"' + '", "'.join(varlist2) + '"\' \n')
    tecplot_macro.append('$!LINEMAP [1]  ASSIGN{YAXISVAR = 11} \n')
    tecplot_macro.append('$!XYLINEAXIS YDETAIL 1 {TITLE{TITLEMODE = USETEXT}} \n')
    tecplot_macro.append(r"$!XYLINEAXIS YDETAIL 1 {TITLE{TEXT = 'Mass Imbalance'}}" + ' \n')
    tecplot_macro.append('$!XYLINEAXIS XDETAIL 1 {TITLE{TITLEMODE = USETEXT}} \n')
    tecplot_macro.append(r"$!XYLINEAXIS XDETAIL 1 {TITLE{TEXT = 'Time'}}" + ' \n')       
    tecplot_macro.append('$!VIEW DATAFIT \n') 
    tecplot_macro.append('$!XYLINEAXIS YDETAIL 1 {RANGEMIN = -1} \n')
    tecplot_macro.append('$!XYLINEAXIS YDETAIL 1 {RANGEMAX = 1} \n')        
    tecplot_macro.append('$!EXPORTSETUP EXPORTFORMAT = PNG \n')
    tecplot_macro.append('$!EXPORTSETUP IMAGEWIDTH = 879 \n')
    tecplot_macro.append('$!EXPORTSETUP EXPORTFNAME = \'' + os.path.join(root, 'convergence.png') + '\' \n')
    tecplot_macro.append('$!EXPORT \n')
    tecplot_macro.append('EXPORTREGION = CURRENTFRAME \n') 
    tecplot_macro.append('$!RemoveVar |MFBD| \n')  
    tecplot_macro.append('$!QUIT \n')
   
    filename = 'ProcessResults.mcr'    
    file_handle = open(filename, 'w')
    file_handle.writelines(tecplot_macro)
    file_handle.close()
   
    print 'Tecplot Macro Written Successfully to ProcessResults.mcr'      
    
if __name__ == "__main__":
    
    # -------------------------
    # --- Default Test Case ---
    # ------------------------- 
    WriteMacro(argv[1])
