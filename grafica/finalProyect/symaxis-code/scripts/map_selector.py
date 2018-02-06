from variable_def import *
from parse_work_directories import *
# Check arguments
if len(sys.argv) < 2 :
    arguments = [ 'python', sys.argv[0], 'work_folder']
    print ' '.join(arguments);
    sys.exit()

print '################################################################'
print '#                         Map Selector                         #'
print '################################################################'

# prepare variables
work_folder = sys.argv[1]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)
results_map = results_folder + '/final.map'

print '<<< Pick the first!'
shutil.copy(results_folder+'/map0/final.map', results_map)
shutil.copy(results_folder+'/map0/tips_align.cor', results_folder+'/tips.cor')
shutil.copy(results_folder+'/axis_align0/align.cor', results_folder+'/axis.cor')
