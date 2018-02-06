from variable_def import *
from parse_work_directories import *
# Check arguments
if len(sys.argv) < 4 :
    arguments = [ 'python', sys.argv[0],'work_folder', 'meshname0', 'meshname1']
    print ' '.join(arguments);
    sys.exit()

def getclassid(meshname):
    candidate = ['mesh', 'cat', 'horse', 'dog', 'centaur', 'david', 'michael', 'victoria', 'wolf', 'gorilla']
    for cl in candidate:
        if meshname.startswith(cl):
            if cl == 'mesh':
                cl = 'SCAPE'
            return cl
    return 'SHREC'    

print '################################################################'
print '#                        Evaluation                            #'
print '################################################################'

###### executables : map2error, TestErrorAllVertices
# prepare variables
work_folder = sys.argv[1]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)
meshnames = [sys.argv[2], sys.argv[3]]
results_map = results_folder + '/final.map'
results_error = results_folder + '/error.txt' 
work_folder = results_folder + '/evaluation'

meshes = [results_folder+'/mesh0.off', results_folder+'/mesh1.off']

# Load configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

cl = [getclassid(meshnames[0]), getclassid(meshnames[1])]
if cl[0] == cl[1] and cl[0] != 'SHREC':
    print 'Check the ground truth of vertex-to-vertex correspondences ...'
    # vertex-to-vertex
    arguments = ['%s/TestErrorAllVertices'%current_path, 
            meshes[0],meshes[1], 
            results_map, 
            work_folder, 
            'default',
            '%s/SurfaceVis'%(external_bin_path),
            bench_settings_path,
            bench_vis_settings_path
            ]
    os.system(' '.join(arguments))
    shutil.copy(work_folder+'/BenchResult_default_mesh0_to_mesh1.txt', results_error)
else:
    # coarse set
    print 'Check the ground truth of sparse correspondences ...'
    benchmark_folder = marked_corrs_path
    arguments = ['%s/map2error'%internal_bin_path, meshes[0], meshes[1], results_map, 
                 benchmark_folder+'/'+meshnames[0]+'.pid', 
                 benchmark_folder+'/'+meshnames[1]+'.pid', results_error, '-geodesic']
    os.system(' '.join(arguments))
