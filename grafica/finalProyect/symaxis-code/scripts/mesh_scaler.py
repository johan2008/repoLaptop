from variable_def import *
######## check and parse arguments
if len(sys.argv) < 4:
    arguments = ['Usage : ', 'python', sys.argv[0], 'fullname-mesh', 'final_mesh', 'work_folder']
    print ' '.join(arguments)
    sys.exit()

print '#################################################################'
print '#                     initialize mesh                           #'
print '#################################################################'

######## executables: msh2msh
full_meshname = sys.argv[1]
final_mesh = sys.argv[2]
work_folder = sys.argv[3]
##### prepare folder
os.system('mkdir -p '+work_folder)

##### scale the mesh to have unit area
arguments = ['%s/msh2msh'%internal_bin_path, full_meshname, final_mesh, '-scale_by_area']
os.system(' '.join(arguments))
