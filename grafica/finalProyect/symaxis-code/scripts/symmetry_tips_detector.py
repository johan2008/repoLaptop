from variable_def import *
from parse_work_directories import *
######## check and parse arguments
if len(sys.argv) < 4:
    arguments = ['Usage : ', 'python', sys.argv[0], 'work_folder', 'index-0-or-1', 'meshname']
    print ' '.join(arguments)
    sys.exit()

print '#################################################################'
print '#                 symmetry and tips detector                    #'
print '#################################################################'

# directory info and meshname
work_folder = sys.argv[1]
index = sys.argv[2]
meshname = sys.argv[3]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)

###### executables : axisrefine, axis2dist, prp2persistence
model = results_folder+'/mesh'+index+'.off'
axis_folder = results_folder+'/axis'+index
symmetry_work_folder = symmetry_work_folder[int(index)]

os.system('mkdir -p '+axis_folder)
os.system('mkdir -p '+symmetry_work_folder)

####### read configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

# check mark mode and copy axes
if config.get('Switches', 'marked_mode') == 'nothing':
    print 'We detect symmetry axes automatically ...'
    arguments = ['python', '%s/symmetry_detector.py'%current_path, work_folder, index]
    os.system(' '.join(arguments))
elif config.get('Switches', 'marked_mode') == 'axis' or config.get('Switches', 'marked_mode') == 'axis_alignment':
    print 'We pick marked symmetry axes!'
    # attempt to copy marked axes
    archive = marked_axis_path + '/'+meshname
    if not os.path.exists(archive):
        print 'We dont find marked axis folder %s, and we have to stop here!'%archive
        sys.exit()
    else:
        os.system('cp -R %s/* %s/'%(archive, axis_folder))

# compute tips
for i in range(3):
    sub_results_folder = axis_folder + '/' + str(i)
    if not os.path.exists(sub_results_folder):
        break
    sub_work_folder = symmetry_work_folder + '/' + str(i)
    os.system('mkdir -p '+sub_work_folder)
    # refine axis
    arguments = ['%s/axisrefine'%internal_bin_path, model, sub_results_folder+'/symmetry.pid', sub_results_folder+'/symmetry.pid']
    os.system(' '.join(arguments))
    # compute distance
    arguments = ['%s/axis2dist'%internal_bin_path, model, sub_results_folder+'/symmetry.pid', 
                 sub_work_folder+'/dist.val', '-geodesic']
    os.system(' '.join(arguments))
    # compute persistent tips
    arguments = ['%s/prp2persistence'%internal_bin_path, model, sub_work_folder+'/dist.val', 
                 '-points', sub_results_folder+'/tips.pid',
                 '-lambda', config.get('TipsDetection-Parameters', 'persistence_lambda')]
    os.system(' '.join(arguments))
