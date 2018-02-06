from variable_def import *
from parse_work_directories import *
######## check and parse arguments
if len(sys.argv) < 4:
    arguments = ['Usage : ', 'python', sys.argv[0], 'fullname-mesh0', 'fullname-mesh1', ' work-folder']
    print ' '.join(arguments)
    sys.exit()
full_name_meshes = [sys.argv[1],sys.argv[2]]
work_folder = sys.argv[3]

####### read configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

####### prepare for variables
meshnames = []
for i in range(2):
    meshname = os.path.splitext(full_name_meshes[i])[0]
    meshname = os.path.basename(meshname)
    meshnames.append(meshname)
pairname = '_'.join(meshnames)


####### prepare for work_folder/results_folder
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)
## mesh related
#os.system('mkdir -p '+mesh_work_folder)
## symmetry related
#os.system('mkdir -p '+symmetry_work_folder)
## axis alignment related
#os.system('mkdir -p '+axisalignment_work_folder)
## mapping related
#os.system('mkdir -p '+map_work_folder)
## evaluation related
#os.system('mkdir -p '+evaluation_work_folder)
# final results
os.system('mkdir -p '+results_folder)
print 'We are runing %s vs %s, working in %s'%(meshnames[0], meshnames[1], work_folder)


####### initialize meshes
if config.get('Pipeline-switches', 'mesh_initializer') == 'on':
    for i in range(2):
        arguments = ['python', '%s/mesh_initializer.py'%current_path, work_folder, str(i), full_name_meshes[i]]
        os.system(' '.join(arguments))

else:
    print 'mesh_initializer is skipped!'

#######
# we have mesh in results_folder/mesh0/mesh.off now
#######


####### compute symmetry
if config.get('Pipeline-switches', 'symmetry_tips_detector') == 'on':
    for  i in range(2):
        arguments = ['python', '%s/symmetry_tips_detector.py'%current_path, work_folder, str(i), meshnames[i]]
        os.system(' '.join(arguments))
else:
    print 'symmetry_tips_detector is skipped!'

#######
# we have everything in results_folder/axis0 now
#######



###### compute axis alignment
if config.get('Pipeline-switches', 'axis_aligner') == 'on':
    arguments = ['python', '%s/axis_aligner.py'%current_path,
                 work_folder, meshnames[0], meshnames[1]]
    os.system(' '.join(arguments))
else:
    print 'axis_aligner is skipped!'

#######
# we have everything in results_folder/axisalign0 and axisalign1
#######



###### compute dense map
if config.get('Pipeline-switches', 'extrapolator') == 'on':
    # here we changed, only extrapolate the first map (with lowest alignment cost)
    for i in range(1):
        arguments = ['python', '%s/extrapolator.py'%current_path,
                     work_folder, str(i),
                     meshnames[0], meshnames[1]]

        os.system(' '.join(arguments))
else:
    print 'extrapolator is skipped!'

######
# we have everything in results_folder/map0 and results_folder/map1
######



###### pick the best map
if config.get('Pipeline-switches', 'map_selector') == 'on':
    # here we changed, simply pick the only map
    arguments = ['python', '%s/map_selector.py'%current_path, work_folder] 
    os.system(' '.join(arguments))
else:
    print 'map_selector is skipped!'

###### evaluate the map
if config.get('Pipeline-switches', 'map_evaluator') == 'on':
    arguments = ['python', '%s/map_evaluator.py'%current_path, work_folder, meshnames[0], meshnames[1]]
    os.system(' '.join(arguments))
else:
    print 'map_evaluator is skipped!'
