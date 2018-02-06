from variable_def import *
from parse_work_directories import *
# Check arguments
if (len(sys.argv) != 4):
    arguments = [ 'python', sys.argv[0] , 'work_folder', 'meshname0', 'meshname1']
    print ' '.join(arguments);
    sys.exit()

print '################################################################'
print '#                       Axis Alignments                        #'
print '################################################################'

###### executables: anchors_to_canonic_align, refinecoarse, axis2prp, myprp2prp, prp_norm, axis2align2, axislength
work_folder = sys.argv[1]
meshnames = [sys.argv[2], sys.argv[3]]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)

def reverse_file(src, dst):
    infile = open(src, 'r')
    lines = infile.readlines()
    infile.close()
    lines.reverse()
    outfile = open(dst, 'w')
    outfile.writelines(lines)
    outfile.close()

# prepare variables
axis_results_folder = [results_folder+'/axis0', results_folder+'/axis1']
align_results_folder = [results_folder+'/axis_align0', results_folder+'/axis_align1']
work_folder = axisalignment_work_folder

# prepare folders
os.system('mkdir -p '+align_results_folder[0])
os.system('mkdir -p '+align_results_folder[1])
os.system('mkdir -p '+work_folder)


# Load configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

# Properties for alignment
properties = [mesh_work_folder[0]+'/cumulative_hist.arff', mesh_work_folder[1]+'/cumulative_hist.arff']
weights = config.get('AxisAlignment-Parameters','curvatures_histogram_weights')

# Paramters
missing_penalty = config.get('AxisAlignment-Parameters', 'missing_penalty')
num_samples = config.get('AxisAlignment-Parameters', 'num_samples')
distortion_thresh = float(config.get('AxisAlignment-Parameters', 'distortion_thresh'))

# Statistics
indices_x = []
indices_y = []
# distortion of the map
distortion_x = []
distortion_y = []
length_axis_x = []
length_axis_y = []
list_length = []
list_cost = []
list_avgcost_inv = []

num_of_axes = [0, 0]
meshes = [results_folder+'/mesh0.off', results_folder+'/mesh1.off']

if config.get('Switches', 'marked_mode') == 'axis_alignment':
    print 'we are generating alignment from marked anchors !'
    alignment_truth_folder = marked_axisalign_path
    anchors = [alignment_truth_folder+'/'+meshnames[0]+'.pid', alignment_truth_folder+'/'+meshnames[1]+'.pid']    
    if not os.path.exists(anchors[0]):
        print 'Unable to find anchors for '+meshname[0]
        sys.exit()
    if not os.path.exists(anchors[1]):
        print 'Unable to find anchors for '+meshname[1]
        sys.exit()
    sub_axis_results_folder0 = '%s/0' % (axis_results_folder[0])
    sub_axis_results_folder1 = '%s/0' % (axis_results_folder[1])

    # generate gt axis alignment
    arguments = ['%s/anchors_to_canonic_align'%internal_bin_path, meshes[0],meshes[1],
                 sub_axis_results_folder0+'/symmetry.pid',sub_axis_results_folder1+'/symmetry.pid',
                 anchors[0], anchors[1], align_results_folder[0]+'/align.cor', '-v']
    os.system(' '.join(arguments))
    # refine the coarse
    arguments = ['%s/refinecoarse'%internal_bin_path, align_results_folder[0]+'/align.cor', align_results_folder[0]+'/align.cor']
    os.system(' '.join(arguments))
    # copy axes
    shutil.copy(sub_axis_results_folder0+'/symmetry.pid', align_results_folder[0]+'/axis0.pid')
    shutil.copy(sub_axis_results_folder1+'/symmetry.pid', align_results_folder[0]+'/axis1.pid')
    shutil.copy(sub_axis_results_folder0+'/symmetry.pid', align_results_folder[1]+'/axis0.pid')
    shutil.copy(sub_axis_results_folder1+'/symmetry.pid', align_results_folder[1]+'/axis1.pid')

    # copy alignments
    shutil.copy(align_results_folder[0]+'/align.cor', align_results_folder[1])

    # copy tips
    shutil.copy(sub_axis_results_folder0+'/tips.pid', align_results_folder[0]+'/tips0.pid')
    shutil.copy(sub_axis_results_folder1+'/tips.pid', align_results_folder[0]+'/tips1.pid')
    shutil.copy(sub_axis_results_folder0+'/tips.pid', align_results_folder[1]+'/tips0.pid')
    shutil.copy(sub_axis_results_folder1+'/tips.pid', align_results_folder[1]+'/tips1.pid')
    
    sys.exit()


def align_in_one_orientation(folder):
    for k in range(2):
        # properties: curvatures (clamped)
        arguments = ['%s/axis2prp'%internal_bin_path, meshes[k], folder+'/axis%d.pid'%k,
                     folder+'/curvatures%d.arff'%k, 
                     '-along 10 0.05 0.10 -across 10 0.05 0.2', '-v']
        os.system(' '.join(arguments))
        arguments = ['%s/myprp2prp'%internal_bin_path, meshes[k],
                     folder+'/curvatures%d.arff'%k, 
                     folder+'/curvatures%d.arff'%k,
                     '-clamp -10 10']
        os.system(' '.join(arguments))
        # properties: geodesic histogram (normalized)
        # normalize geodesic
        arguments = ['%s/prp_norm'%internal_bin_path, meshes[k], 
                     properties[k],
                     folder+'/geodesic%d.arff'%k]
        os.system(' '.join(arguments))
        # combine geodesic and curvatures
        arguments = ['%s/myprp2prp'%internal_bin_path, meshes[k], 
                     folder+'/geodesic%d.arff'%k,
                     folder+'/final%d.arff'%k,
                     '-add ', folder+'/curvatures%d.arff'%k, weights, '-v']
        os.system(' '.join(arguments))

    arguments = ['%s/axis2align2'%internal_bin_path,
                 meshes[0], meshes[1],
                 folder+'/axis0.pid', 
                 folder+'/axis1.pid',
                 '-property', folder+'/final0.arff', 
                 folder+'/final1.arff',
                 '%s/align.cor'% (folder),
                 '-start', '1',
                 '-missing', missing_penalty,
                 '-score_file', '%s/score.txt'%(folder), 
                 '-nsamples', num_samples,
                 '-length_file', '%s/length.txt'%(folder)]
    os.system(' '.join(arguments))


# for all combination
# we assume there are at most 3 axes for each mesh
for i in range(5):
    # Check the existence of axis subdir
    sub_axis_results_folder0 = '%s/%d' % (axis_results_folder[0], i)
    if not os.path.exists(sub_axis_results_folder0):
        num_of_axes[0] = i
        break
    for j in range(5):
        # Check the existence of axis subdir
        sub_axis_results_folder1 = '%s/%d' % (axis_results_folder[1], j)
        if not os.path.exists(sub_axis_results_folder1):
            num_of_axes[1] = j
            break
        print 'Aligning %d_%d ...'%(i,j)
        sub_work_folder = ['%s/%d_%d/head' % (work_folder, i, j),
                              '%s/%d_%d/tail' % (work_folder, i, j)]

        # original orientation
        commands = ['mkdir', '-p', sub_work_folder[0]]
        os.system(' '.join(commands))
        # copy axes
        shutil.copy(sub_axis_results_folder0+'/symmetry.pid', sub_work_folder[0]+'/axis0.pid')
        shutil.copy(sub_axis_results_folder1+'/symmetry.pid', sub_work_folder[0]+'/axis1.pid')
        align_in_one_orientation(sub_work_folder[0])
        # reversed orientation
        commands = ['mkdir', '-p', sub_work_folder[1]]
        os.system(' '.join(commands))
        # copy axes (and reverse the second)
        shutil.copy(sub_axis_results_folder0+'/symmetry.pid', sub_work_folder[1]+'/axis0.pid')
        reverse_file(sub_axis_results_folder1+'/symmetry.pid', sub_work_folder[1]+'/axis1.pid')
        align_in_one_orientation(sub_work_folder[1])

        indices_x.append(str(i))
        indices_y.append(str(j))
        # distortion of the self-maps
        if os.path.exists('%s/distortion.txt'%sub_axis_results_folder0):
            fin = open('%s/distortion.txt'%sub_axis_results_folder0)
            distortion_x.append(float(fin.readline()))
            fin.close()
        else:
            distortion_x.append(1.0)
        if os.path.exists('%s/distortion.txt'%sub_axis_results_folder1):
            fin = open('%s/distortion.txt'%sub_axis_results_folder1)
            distortion_y.append(float(fin.readline()))
            fin.close()
        else:
            distortion_y.append(1.0)

        fin = open(sub_work_folder[0]+'/score.txt')
        score0 = float(fin.readline())
        if score0 < 0.0001:
            score0 = 0.0001
        fin.close()
        fin = open(sub_work_folder[1]+'/score.txt')
        score1 = float(fin.readline())
        if score1 < 0.0001:
            score1 = 0.0001
        fin.close()
        if score0 > score1:
            score = score1
        else:
            score = score0
        fin = open(sub_work_folder[0]+'/length.txt')
        length0 = int(fin.readline())
        fin.close()
        fin = open(sub_work_folder[1]+'/length.txt')
        length1 = int(fin.readline())
        fin.close()
        if length0 > length1:
            length = length0
        else:
            length = length1
        list_length.append(length)
        list_cost.append(score)
        avgcost0_inv = length0/score0
        avgcost1_inv = length1/score1
        if avgcost0_inv > avgcost1_inv:
            avgcost_inv = avgcost0_inv
        else:
            avgcost_inv = avgcost1_inv
        list_avgcost_inv.append(avgcost_inv)

        # Add axis length
        fin = open(sub_work_folder[0]+'/axis0.pid')
        lines = fin.read()
        lines = lines.split('\n')
        length_axis_x.append(len(lines)-1)
        fin.close()
        fin = open(sub_work_folder[0]+'/axis1.pid')
        lines = fin.read()
        lines = lines.split('\n')
        length_axis_y.append(len(lines)-1)
        fin.close()
max_length = 0
min_cost = 1000000.0
max_sym_sym_alignlen_cost = 0
max_sym_sym_cost = 0
max_sym_sym_axislen_cost = 0
best_x = -1
best_y = -1

pick_align_method = config.get('AxisAlignment-Parameters', 'best_alignment_strategy')

# whether we check the length of axes?
if num_of_axes[0] > 1 and num_of_axes[1] > 1:
    axis_length_checking_flag = 1
else:
    print 'One mesh has only one symmetry axis, and we will not check axis length'
    axis_length_checking_flag = 0
length_x = []
length_y = []
length_thresh = config.get('AxisAlignment-Parameters', 'axis_length_thresh')
if axis_length_checking_flag == 1:
    for i in range(num_of_axes[0]):
        arguments = ['%s/axislength'%internal_bin_path, meshes[0], 
                     '%s/%d/symmetry.pid' % (axis_results_folder[0], i), '%s/axis_length_temp.txt'%work_folder]
        os.system(' '.join(arguments))
        infile = open('%s/axis_length_temp.txt'%work_folder, 'r')
        line = infile.read()
        infile.close()
        length_x.append(float(line))
    for i in range(num_of_axes[1]):
        arguments = ['%s/axislength'%internal_bin_path, meshes[1], 
                     '%s/%d/symmetry.pid'%(axis_results_folder[1], i), '%s/axis_length_temp.txt'%work_folder]
        os.system(' '.join(arguments))
        infile = open('%s/axis_length_temp.txt'%work_folder, 'r')
        line = infile.read()
        infile.close()
        length_y.append(float(line))
else:
    for i in range(num_of_axes[0]):
        length_x.append(10.0)
    for i in range(num_of_axes[1]):
        length_y.append(10.0)

print 'Summary :'
print 'indices_x = '+str(len(indices_x))
print 'indices_y = '+str(len(indices_y))
print 'length_x = '+str(len(length_x))
print 'length_y = '+str(len(length_y))



for i in range(len(indices_x)):
    eval_sym_sym_alignlen_cost = (float(distortion_x[i])-distortion_thresh)*(float(distortion_y[i])-distortion_thresh)*float(list_avgcost_inv[i])
    eval_sym_sym_axislen_cost = (float(distortion_x[i])-distortion_thresh)*(float(distortion_y[i])-distortion_thresh)*float(length_axis_x[i])*float(length_axis_y[i])/float(list_cost[i])
    eval_sym_sym_cost =  (float(distortion_x[i])-distortion_thresh)*(float(distortion_y[i])-distortion_thresh)/float(list_cost[i])
#    eval_sym_sym_cost =  float(distortion_x[i])*float(distortion_y[i])/float(list_cost[i])
    print '%d: axis(%d,%d) distortion(%f, %f) length_of_axes(%d %d) length_of_alignment(%d) cost(%f)'%(i,int(indices_x[i]),int(indices_y[i]),float(distortion_x[i]),float(distortion_y[i]), int(length_axis_x[i]), int(length_axis_y[i]), int(list_length[i]),float(list_cost[i]))
    print '        eval_sym_sym_alignlen_cost(%f)'%(eval_sym_sym_alignlen_cost)
    print '        eval_sym_sym_axislen_cost(%f)'%(eval_sym_sym_axislen_cost)
    print '        eval_sym_sym_cost(%f)'%(eval_sym_sym_cost)
    if length_x[int(indices_x[i])] < float(length_thresh) or length_y[int(indices_y[i])] < float(length_thresh):
        print '     axis length is shorter than thresholds'
        continue
    if pick_align_method == 'lowest-cost' and list_cost[i]<min_cost:
        best_x = indices_x[i]
        best_y = indices_y[i]
        min_cost = list_cost[i]
    elif pick_align_method == 'longest' and list_length[i]>max_length:
        best_x = indices_x[i]
        best_y = indices_y[i]
        max_length = list_length[i]
    elif pick_align_method == 'sym_sym_alignlen_cost' and eval_sym_sym_alignlen_cost > max_sym_sym_alignlen_cost:
        best_x = indices_x[i]
        best_y = indices_y[i]
        max_sym_sym_alignlen_cost = eval_sym_sym_alignlen_cost
    elif pick_align_method == 'sym_sym_axislen_cost' and eval_sym_sym_axislen_cost > max_sym_sym_axislen_cost:
        best_x = indices_x[i]
        best_y = indices_y[i]
        max_sym_sym_axislen_cost = eval_sym_sym_axislen_cost
    elif pick_align_method == 'sym_sym_cost' and eval_sym_sym_cost > max_sym_sym_cost:
        best_x = indices_x[i]
        best_y = indices_y[i]
        max_sym_sym_cost = eval_sym_sym_cost
print 'Best alignment strategy: %s'%(pick_align_method)
print 'Best alignment: %s %s'%(best_x, best_y)

best_x = int(best_x)
best_y = int(best_y)
# copy axes
shutil.copy('%s/%d_%d/head/axis0.pid' % (work_folder, best_x, best_y), align_results_folder[0])
shutil.copy('%s/%d_%d/head/axis1.pid' % (work_folder, best_x, best_y), align_results_folder[0])
shutil.copy('%s/%d_%d/tail/axis0.pid' % (work_folder, best_x, best_y), align_results_folder[1])
shutil.copy('%s/%d_%d/tail/axis1.pid' % (work_folder, best_x, best_y), align_results_folder[1])

# copy alignments
shutil.copy('%s/%d_%d/head/align.cor' % (work_folder, best_x, best_y), align_results_folder[0])
shutil.copy('%s/%d_%d/tail/align.cor' % (work_folder, best_x, best_y), align_results_folder[1])

# copy scores
shutil.copy('%s/%d_%d/head/score.txt' % (work_folder, best_x, best_y), align_results_folder[0])
shutil.copy('%s/%d_%d/tail/score.txt' % (work_folder, best_x, best_y), align_results_folder[1])

# copy tips
shutil.copy('%s/%d/tips.pid'%(axis_results_folder[0], best_x), 
            align_results_folder[0]+'/tips0.pid')
shutil.copy('%s/%d/tips.pid'%(axis_results_folder[1], best_y), 
            align_results_folder[0]+'/tips1.pid')
shutil.copy('%s/%d/tips.pid'%(axis_results_folder[0], best_x), 
            align_results_folder[1]+'/tips0.pid')
shutil.copy('%s/%d/tips.pid'%(axis_results_folder[1], best_y), 
            align_results_folder[1]+'/tips1.pid')

# copy symmetry.map
if os.path.exists('%s/%d/symmetry.map'%(axis_results_folder[0], best_x)):
    shutil.copy('%s/%d/symmetry.map'%(axis_results_folder[0], best_x), align_results_folder[0]+'/symmetry0.map')
if os.path.exists('%s/%d/symmetry.map'%(axis_results_folder[1], best_y)):
    shutil.copy('%s/%d/symmetry.map'%(axis_results_folder[1], best_y), align_results_folder[0]+'/symmetry1.map')
if os.path.exists('%s/%d/symmetry.map'%(axis_results_folder[0], best_x)):
    shutil.copy('%s/%d/symmetry.map'%(axis_results_folder[0], best_x), align_results_folder[1]+'/symmetry0.map')
if os.path.exists('%s/%d/symmetry.map'%(axis_results_folder[1], best_y)):
    shutil.copy('%s/%d/symmetry.map'%(axis_results_folder[1], best_y), align_results_folder[1]+'/symmetry1.map')

# swap align_results_folder if necessary
scores = []
infile = open(align_results_folder[0]+'/score.txt', 'r')
val = float(infile.read())
infile.close()
scores.append(val)

infile = open(align_results_folder[1]+'/score.txt', 'r')
val = float(infile.read())
infile.close()
scores.append(val)
print '*****Final competition : %.3f %.3f'%(scores[0], scores[1])
if scores[0] > scores[1]:
    # swap two folders
    print 'Swap the folders!'
    os.system('mv '+align_results_folder[0]+' '+align_results_folder[0]+'_temp')
    os.system('mv '+align_results_folder[1]+' '+align_results_folder[0])
    os.system('mv '+align_results_folder[0]+'_temp '+align_results_folder[1])
