from variable_def import *
from parse_work_directories import *
######## check and parse arguments
if len(sys.argv) < 4:
    arguments = ['Usage : ', 'python', sys.argv[0], 'work_folder', 'index-0-or-1', 'full-mesh-name']
    print ' '.join(arguments)
    sys.exit()

print '#################################################################'
print '#                     initialize mesh                           #'
print '#################################################################'

######## executables: msh2msh
work_folder = sys.argv[1]
index = sys.argv[2]
full_meshname = sys.argv[3]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)
final_mesh = results_folder+'/mesh'+index+'.off'
mesh_work_folder = mesh_work_folder[int(index)]
##### prepare folder
os.system('mkdir -p '+mesh_work_folder)

# Load configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

##### scale the mesh to have unit area
arguments = ['%s/msh2msh'%internal_bin_path, full_meshname, final_mesh, '-scale_by_area']
os.system(' '.join(arguments))

##### compute histogram
num_samples = int(config.get('AxisAlignment-Parameters', 'num_histogram_samples'))
arguments = ['%s/msh2prp'%internal_bin_path, full_meshname, mesh_work_folder+'/cumulative_hist.arff', '-dijkstra_histogram', '-v', '-num_samples %d'%num_samples]
os.system(' '.join(arguments))
