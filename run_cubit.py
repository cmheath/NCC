import cubit
from sys import argv
import subprocess
import numpy as np

def check_quality(volume):

    proc = subprocess.Popen(['tail', '-10', 'python.stdout'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    out, err = proc.communicate() 
    for line in out.split('\n'): 
        if '1 volume(s) did not mesh' in line:
            return -1

    cubit.cmd('quality volume ' + str(volume) + ' Aspect Ratio global')
    proc = subprocess.Popen(['tail', '-10', 'python.stdout'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    out, err = proc.communicate() 
    for line in out.split('\n'): 
        if 'Aspect' not in line: 
            continue 
        fields = line.split()
        print 'Max Aspect Ratio = ' + fields[-2]
    return float(fields[-2])

def find_closest_surfaces(centroids_t):

    num_surfaces = cubit.get_surface_count()
    centroids = []
    
    for i in range(1, num_surfaces+1):
        centroids.append(np.array(cubit.get_surface_centroid(i)))
    
    dist = [[0 for x in xrange(len(centroids))] for x in xrange(len(centroids_t))]
    
    for i, pt1 in enumerate(centroids):
        for j, pt2 in enumerate(centroids_t):
            dist[j][i] = np.linalg.norm(pt2 - pt1)

    dist =  np.array(np.transpose(dist))
    
    surface_ids = ''
    
    for id in dist.argmin(axis=0):
        surface_ids += (str(id+1) + ' ')
        
    return (surface_ids)
        
dome_plane = float(argv[2])
venturi_angle = float(argv[3])
vane_height = float(argv[4])

cubit.init(" ")

# --- Create New Cubit File
cubit.cmd('new')

# --- Import STEP file with automatic heal
fname = 'Config' + argv[1] + '.step'
cubit.cmd('import step "%s" heal'%fname)

cubit.cmd('Set Max Memory On 500') 

cubit.cmd('regularize volume all')
cubit.cmd('undo on')    
cubit.cmd('set multisweep off') 

p0 = -3.0
p1 = -2.5
p2 = -2.0 # --- Do not edit this line (corresponds to exact geometry)
p3 = dome_plane + 0.2
p4 = p3 + 0.4
p5 = p4 + 0.5
p6 = p5 + 1.0
p7 = 7.5

cubit.cmd('webcut volume all with plane xplane offset ' + str(dome_plane) + ' imprint merge ')

cubit.cmd('compress ids all')

cubit.cmd('create pressure 1 on surface ( at ' + str(p0) + ' 9.229457e+000 2.032004e-001 ordinal 1 ordered) magnitude 1')
cubit.cmd('create pressure 2 on surface ( at ' + str(p7) + ' 9.443297e+000 2.761321e-001 ordinal 1 ordered) magnitude 2')

centroids_t = []
centroids_t.append(np.array([(3.515369 + dome_plane / 2.0), 10.74310, 1.129146]))
centroids_t.append(np.array([(3.7505 + dome_plane / 2.0), 10.04726, 1.556210]))
centroids_t.append(np.array([(3.7505 + dome_plane / 2.0), 8.453696, 1.388950]))
centroids_t.append(np.array([(3.093733 + dome_plane / 2.0), 7.709811, 0.8103338]))
centroids_t.append(np.array([3.75, 9.253586, 0.9725911]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 10.74310, 1.129146]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 10.04614, 1.556205]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 9.253586, 0.9725911]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 8.454910, 1.388959]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 7.709811, 0.8103338]))
centroids_t.append(np.array([((p3 - dome_plane) / 2.0 + dome_plane), 9.253586, 0.9725911]))

surface_ids = find_closest_surfaces(np.array(centroids_t))
cubit.cmd('create pressure 3 on surface ' + surface_ids + 'magnitude 3')

centroids_t = []
centroids_t.append(np.array([(3.515369 + dome_plane / 2.0), 10.74310, -1.129146]))
centroids_t.append(np.array([3.750500, 10.05354, -0.5562327]))
centroids_t.append(np.array([3.750000, 9.253586, -0.9725911]))
centroids_t.append(np.array([3.750500, 8.459906, -0.3889721]))
centroids_t.append(np.array([(3.093733 + dome_plane / 2.0), 7.709811, -0.8103338]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 10.74310, -1.129146]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 10.05226, -0.5562232]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 9.253586, -0.9725911]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 8.461027, -0.3889777]))
centroids_t.append(np.array([((p2 - p0) / 2.0 + p0), 7.709811, -0.8103338]))
centroids_t.append(np.array([((p3 - dome_plane) / 2.0 + dome_plane), 9.253586, -0.9725911]))

surface_ids = find_closest_surfaces(np.array(centroids_t))
cubit.cmd('create pressure 4 on surface ' + surface_ids + 'magnitude 4')
                                               
centroids_t = []
centroids_t.append(np.array([4.688935, 10.91725, -0.5859262]))
centroids_t.append(np.array([4.689590, 10.90707, 0.7512378]))
centroids_t.append(np.array([5.810667, 7.795583, -0.4340362]))
centroids_t.append(np.array([5.805872, 7.775533, 0.6772361]))
centroids_t.append(np.array([(0.9875 + dome_plane / 2.0), 11.00000, 0.0]))
centroids_t.append(np.array([(2.020312 + dome_plane / 2.0), 7.50000, 0.0]))

surface_ids = find_closest_surfaces(np.array(centroids_t))
cubit.cmd('create pressure 5 on surface ' + surface_ids + 'magnitude 5')
 
# --- Subdivide volume    
cubit.cmd('webcut volume all with plane xplane offset ' + str(p1) + ' imprint merge ')
cubit.cmd('webcut volume all with plane xplane offset ' + str(p2) + ' imprint merge ')

pilot_cut = 0.0 - ((1.15/2 - 0.55/2)/np.tan(np.pi*venturi_angle/180.0) + (0.55/2 - 0.03/2)/np.tan(np.pi*venturi_angle/180.0) + vane_height + 0.03 + 0.15)
stage2_cut = dome_plane - ((0.85/2 - 0.55/2)/np.tan(np.pi*venturi_angle/180.0) + (0.55/2 - 0.03/2)/np.tan(np.pi*venturi_angle/180.0) + vane_height + 0.03 + 0.15)
stage3_cut = dome_plane - ((0.85/2 - 0.55/2)/np.tan(np.pi*venturi_angle/180.0) + (0.55/2 - 0.03/2)/np.tan(np.pi*venturi_angle/180.0) + vane_height + 0.03 + 0.15)
cubit.cmd('create curve arc radius .5 center location ' + str(stage2_cut - 0.25) + ' 10.425 0.000000 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id))
cubit.cmd('create curve arc radius .5 center location ' + str(stage2_cut - 0.25) + ' 10.0492037 1.05621387 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id))  
cubit.cmd('create curve arc radius .5 center location ' + str(stage3_cut - 0.25) + ' 8.00000 0.000000 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id)) 
cubit.cmd('create curve arc radius .5 center location ' + str(stage3_cut - 0.25) + ' 8.457968 0.8889683 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id))
cubit.cmd('create curve arc radius .5 center location ' + str(stage2_cut) + ' 10.425 0.000000 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id))
cubit.cmd('create curve arc radius .5 center location ' + str(stage2_cut) + ' 10.0492037 1.05621387 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id))  
cubit.cmd('create curve arc radius .5 center location ' + str(stage3_cut) + ' 8.00000 0.000000 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id)) 
cubit.cmd('create curve arc radius .5 center location ' + str(stage3_cut) + ' 8.457968 0.8889683 normal 1 0 0 start angle 0 stop angle 360')
last_id = cubit.get_last_id("curve")
cubit.cmd('webcut volume all with loop curve ' + str(last_id) + ' imprint merge')
cubit.cmd('delete curve ' + str(last_id))

cubit.cmd('webcut volume all with plane xplane offset ' + str(p3) + ' imprint merge')
cubit.cmd('webcut volume all with plane xplane offset ' + str(p4) + ' imprint merge')
cubit.cmd('webcut volume all with plane xplane offset ' + str(p5) + ' imprint merge')
cubit.cmd('webcut volume all with plane xplane offset ' + str(p6) + ' imprint merge')                                               

# --- Mesh     
cubit.cmd('volume all scheme tetmesh')

scale = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 7]
min_depth = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
max_depth = [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
min_num_layers_3D = [4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1]
# max_size_factor = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 2.25, 2.25, 2.25, 2.25, 2.25, 3.125, 3.125, 4.5, 6]
max_size_factor = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 2.25, 2.25, 2.25, 2.25, 1.75, 3.125, 2.75, 4.5, 6]
volume = [8, 7, 6, 5, 4, 17, 13, 14, 15, 16, 18, 9, 10, 11, 12, 2, 19, 3, 20, 1]

scale_factor = 1.025

for index in range(20):
    max_size = max_size_factor[index] * 0.0315
    min_size = max_size * 0.3
    max_gradient = 1.15
    min_ar = 100 # --- Initialize min_ar
    
    cubit.cmd('volume ' + str(volume[index]) + ' sizing function type skeleton scale ' + str(scale[index]) + ' time_accuracy_level 2.0 min_size ' + str(min_size) + ' max_size ' + str(max_size) + ' max_gradient ' + str(max_gradient) + \
              ' min_depth ' + str(min_depth[index]) + ' max_depth ' + str(max_depth[index]) + ' min_num_layers_3d ' + str(min_num_layers_3D[index]) + ' min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')         
    cubit.cmd('mesh volume ' + str(volume[index]))
    
    ar = check_quality(volume[index])
    if ar < min_ar and ar != -1:
        min_ar = ar
        max_size_stored = max_size
        max_gradient_stored = max_gradient
        scale_stored = scale[index]
        print "Updating best case: max size: %f , max_grad: %f , scale: %f" % (max_size_stored, max_gradient_stored, scale_stored) 
        
    if ar > 6.0 and ar < 50.0: # --- Attempt condition number smoothing
        cubit.cmd('volume ' + str(volume[index]) + ' smooth scheme condition number beta 1.5 cpu 2.0')    
        cubit.cmd('smooth volume ' + str(volume[index]))

    for val in range(3):
        scale[index] *= scale_factor
        
        ar = check_quality(volume[index])
        if ar < min_ar and ar != -1:
            min_ar = ar
            max_size_stored = max_size
            max_gradient_stored = max_gradient
            scale_stored = scale[index]
            print "Updating best case: max size: %f , max_grad: %f , scale: %f" % (max_size_stored, max_gradient_stored, scale_stored)            
        
        while (ar > 6.0 and max_gradient < 2.0) or ar == -1:
            cubit.cmd('delete mesh volume ' + str(volume[index]) + ' propagate')
            max_gradient += max_gradient * 0.025

            cubit.cmd('volume ' + str(volume[index]) + ' sizing function type skeleton scale ' + str(scale[index]) + ' time_accuracy_level 2.0 min_size ' + str(min_size) + ' max_size ' + str(max_size) + ' max_gradient ' + str(max_gradient) + \
                      ' min_depth ' + str(min_depth[index]) + ' max_depth ' + str(max_depth[index]) + ' min_num_layers_3d ' + str(min_num_layers_3D[index]) + ' min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')         
            cubit.cmd('mesh volume ' + str(volume[index]))
            
            ar = check_quality(volume[index]) 
            if ar < min_ar and ar != -1:
                min_ar = ar
                max_size_stored = max_size
                max_gradient_stored = max_gradient
                scale_stored = scale[index]
                print "Updating best case: max size: %f , max_grad: %f , scale: %f" % (max_size_stored, max_gradient_stored, scale_stored)                
            
            if ar > 6.0 and ar < 50.0: # --- Attempt condition number smoothing
                cubit.cmd('volume ' + str(volume[index]) + ' smooth scheme condition number beta 1.5 cpu 2.0')    
                cubit.cmd('smooth volume ' + str(volume[index])) 

            ar = check_quality(volume[index])
            if ar < min_ar and ar != -1:
                min_ar = ar
                max_size_stored = max_size
                max_gradient_stored = max_gradient
                scale_stored = scale[index]
                print "Updating best case: max size: %f , max_grad: %f , scale: %f" % (max_size_stored, max_gradient_stored, scale_stored)
        
        if ar < 6.0 and ar != -1:
            break

    if ar > 6.0: # Return to best solution identified prior
        cubit.cmd('delete mesh volume ' + str(volume[index]) + ' propagate')
        
        max_gradient = max_gradient_stored
        max_size = max_size_stored
        min_size = max_size * 0.3
        scale[index] = scale_stored    

        cubit.cmd('volume ' + str(volume[index]) + ' sizing function type skeleton scale ' + str(scale[index]) + ' time_accuracy_level 2.0 min_size ' + str(min_size) + ' max_size ' + str(max_size) + ' max_gradient ' + str(max_gradient) + \
                  ' min_depth ' + str(min_depth[index]) + ' max_depth ' + str(max_depth[index]) + ' min_num_layers_3d ' + str(min_num_layers_3D[index]) + ' min_num_layers_2d 1 min_num_layers_1d 1 max_span_ang_surf 45 max_span_ang_curve 45')         
        cubit.cmd('mesh volume ' + str(volume[index]))
        ar = check_quality(volume[index])            
        if ar > 6.0 and ar < 50.0: # --- Attempt condition number smoothing
            cubit.cmd('volume ' + str(volume[index]) + ' smooth scheme condition number beta 1.5 cpu 2.0')    
            cubit.cmd('smooth volume ' + str(volume[index]))

cubit.cmd('volume all smooth scheme condition number beta 2.0 cpu 4.0')
cubit.cmd('smooth volume all')

# --- Save and export Patran             
cubit.cmd('save as "Config' + argv[1] + '.cub" overwrite')                                                             
cubit.cmd('export patran "patran.out" overwrite')

cubit.cmd('reset')
