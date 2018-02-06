import time,os,sys,ConfigParser,shutil
from parse_work_directories import *
from gen_extrinsic_sym import *
from gen_intrinsic_sym import *
def max_dec(val_x, val_y):
    if val_x > val_y:
        return -1
    else:
        return 1
if len(sys.argv) < 3:
    arguments = [ 'python', sys.argv[0], 'work_folder', 'index-0-or-1']
    print ' '.join(arguments)
    sys.exit()

# directory info and meshname
work_folder = sys.argv[1]
index = sys.argv[2]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)

model = results_folder+'/mesh'+index+'.off'
axis_folder = results_folder+'/axis'+index
symmetry_work_folder = symmetry_work_folder[int(index)]
mesh_work_folder = mesh_work_folder[int(index)]

os.system('mkdir -p '+axis_folder)
os.system('mkdir -p '+symmetry_work_folder)

####### read configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

# generate tips and agds
tips = '%s/sym-tips.pid'%mesh_work_folder
agd = '%s/agd.val'%mesh_work_folder
print 'generating tips and agds ...'
tips_agd_work_folder = mesh_work_folder + '/agd'
os.system('mkdir -p '+tips_agd_work_folder)
arguments = ['%s/IntrinsicSymmetryBlended'%current_path, 
             model, tips_agd_work_folder,
            '%s/SurfaceVis'%(external_bin_path),
             agd_settings_path,
            agd_vis_settings_path]
os.system(' '.join(arguments))
agd_temp = '%s/mesh%s/IS_Blended_mesh%s.feat.vals'%(tips_agd_work_folder,index, index)
tips_temp = '%s/mesh%s/IS_Blended_mesh%s.smplset.agd'%(tips_agd_work_folder, index, index)
arguments = ['%s/vprp2prp'%internal_bin_path, model, agd_temp, agd]
os.system(' '.join(arguments))
arguments = ['%s/vpts2pts'%internal_bin_path, model, tips_temp, tips]
os.system(' '.join(arguments))

# extrinsic
max_extrinsic_error = float(config.get('SymDetection-Parameters', 'extrinsic_max_error'))
start_extrinsic = time.time()
work_dir = symmetry_work_folder + '/extrinsic'
extrinsic_suc_flag = GenExtrinsicSym(model, work_dir, axis_folder, max_extrinsic_error)
elapsed_extrinsic = (time.time() - start_extrinsic)

if extrinsic_suc_flag > 0:
    flag_file = open('%s/flag.txt'%(axis_folder), 'w')
    flag_file.write('extrinsic')
    flag_file.close()
    sys.exit()
else:
    print 'Extrinsic method fails. Let\'s start intrinsic method!'

# intrinsic
start_intrinsic = time.time()
work_dir = symmetry_work_folder + '/intrinsic'
intrinsic_suc_flag = GenIntrinsicSym(model, tips, agd, work_dir, axis_folder, config)
elapsed_intrinsic = (time.time() - start_intrinsic)
print '************Intrinsic total running time (excluding selection): %f seconds'%elapsed_intrinsic
if intrinsic_suc_flag:
    flag_file = open('%s/flag.txt'%(axis_folder), 'w')
    flag_file.write('intrinsic')
    flag_file.close()
    sys.exit()
else:
    print 'Unable to find any valid axis!'

#in_work_dir = 'WorkDir/%s'%(meshname)
#errors_info_file = open('%s/IS_Blended_%sIS_Blended_DenseMapPreceiseAll__%s_to_%s.dense.preceise.confidences.txt'%(in_work_dir, meshname, meshname, meshname), 'r')
#line = errors_info_file.readline()
#num_axes_actual = line.split(' ')
#num_axes_actual = int(num_axes_actual[1])
#num_candidates = int(config.get('Intrinsic', 'CandidateAxes'))
#if num_axes_actual > num_candidates:
    #num_axes_actual = num_candidates
#line = errors_info_file.readline()
#errors_info_file.close()
#error_cands = line.split(' ')
#num_in_axes = int(config.get('Intrinsic', 'NumAxes'))
##if num_in_axes < num_axes_actual:
##    num_axes_actual = num_in_axes
## another i, avoiding zero-length axes
##counter_i = 0
## tips, for checking the consistency
#tips = '%s/%s.pid'%(config.get('Directories', 'TipsFolder'), meshname)
## Keep all legal axes
## evaluation = (distortion - DistortionThresh) * axis_length
#thresh = float(config.get('Intrinsic', 'DistortionThresh'))
#evaluation = []
#indices = []
#errors = []
#axislen = []
#valid_counter = 0
#for i in range(num_axes_actual):
    #error = float(error_cands[2*i])
    #errors.append(error)
    #in_axes = '%s/%d/symmetry.pid'%(in_work_dir, i)
    #in_syms = '%s/%d/symmetry.val'%(in_work_dir, i)

    #in_axes_file = open(in_axes, 'r');
    #lines = in_axes_file.read()
    #num_lines = len(lines.split('\n')) - 1
    #axislen.append(num_lines)

    #indices.append(i)
    ## check whether the axis is empty
    #if num_lines <= 0:
        #evaluation.append(0.0)
    #else:
        ## check whether the axis is consistent to the val
        #consistency = '%s/%d/consistency.txt'%(in_work_dir, i)
        #commands = ['check_axis_consistency', mesh, tips, in_axes, in_syms, consistency,
                #'-thres_axis', '0.33' ]
        #os.system(' '.join(commands))
        #in_consistency_file = open(consistency, 'r')
        #lines = in_consistency_file.readline()
        #if lines.strip() == 'no':
            #evaluation.append(0.0)
        #else:
            #evaluation.append((error - thresh)*float(num_lines))
            #valid_counter = valid_counter + 1
    #print evaluation[len(evaluation)-1]
## if we cannot find any axis, we just add the first unempty axis
#if valid_counter == 0 and num_axes_actual >= 1:
    #print 'add the first one back. (all axes on not consistent!)'
    #evaluation[0] = (errors[0] - thresh)*float(axislen[0])

## sort evaluations
#indices.sort(lambda x,y: max_dec(evaluation[x], evaluation[y]))

#if num_in_axes > len(indices):
    #num_in_axes = len(indices)
#for i in range(num_in_axes):
    #print evaluation[indices[i]]
    #if evaluation[indices[i]] <= 0.0:
        #continue
    #print '     Axis %d: %d (confidence: %f, length: %d, eval: %f)'%(i, indices[i], errors[indices[i]], axislen[indices[i]], evaluation[indices[i]])
    #print '     copying ...'
    #syms_sub_dir = '%s/%d'%(syms_dir, i)
    #in_axes = '%s/%d/symmetry.pid'%(in_work_dir, indices[i])
    #in_maps = '%s/%d/symmetry.map'%(in_work_dir, indices[i])

    #os.system('mkdir -p %s'%(syms_sub_dir))
    #shutil.copy(in_axes, syms_sub_dir)
    #shutil.copy(in_maps, syms_sub_dir)
    #outfile = open('%s/distortion.txt'%(syms_sub_dir), 'w')
    #outfile.write(str(errors[indices[i]]))
    #outfile.close()

#os.system('mkdir -p %s'%(syms_dir))
#flag_file = open('%s/flag.txt'%(syms_dir), 'w')
#flag_file.write('intrinsic')
#flag_file.close()
#sub_archive_dir = '%s/%s'%(archive_dir, meshname)
#if os.path.exists(sub_archive_dir):
    #os.system('rm -f -r %s'%(sub_archive_dir))
#os.system('mkdir -p %s'%(sub_archive_dir))
#os.system('cp -r %s/* %s'%(syms_dir, sub_archive_dir))
