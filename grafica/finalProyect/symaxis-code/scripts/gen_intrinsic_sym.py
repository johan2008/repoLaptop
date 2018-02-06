from variable_def import *
# executables: intrinsic_triplets_selection, precisemap2vertexmap, symmap2dist, symdist2binary, symdist2axis, axisrefine,check_axis_consistency
def GetAxisLength(axis):
    fin = open(axis, 'r')
    lines = fin.readlines()
    fin.close()
    return len(lines)-1

# call by sorting
def max_dec(val_x, val_y):
    if val_x > val_y:
        return -1
    else:
        return 1

def IsAxisValid(mesh, tips, axis, sym, work_dir):
    axis_length = GetAxisLength(axis)
    if axis_length <= 0:
        return False
    consistent_flag = '%s/consistent_flag.txt'%work_dir
    commands = ['%s/check_axis_consistency'%internal_bin_path, mesh, tips, axis, sym, consistent_flag, '-thres_axis', '0.33' ]
    os.system(' '.join(commands))
    fin = open(consistent_flag, 'r')
    line = fin.readline()
    if line.strip() == 'no':
        return False
    return True

def GenIntrinsicSym(mesh, tips, agd, work_dir, result_dir, config):
    # mkdir
    os.system('mkdir -p '+work_dir)
    # parameters
    min_num_shared = config.get('SymDetection-Parameters', 'intrinsic_min_num_shared')
    min_dist_thresh = config.get('SymDetection-Parameters', 'intrinsic_min_dist_thresh')
    dist_to_other_thresh = config.get('SymDetection-Parameters', 'intrinsic_dist_to_other_thresh')
    feature_thresh = config.get('SymDetection-Parameters', 'intrinsic_feature_thresh')
    map_to_dist_max = config.get('SymDetection-Parameters', 'intrinsic_map_to_dist_max')
    num_selected_axes = int(config.get('SymDetection-Parameters', 'intrinsic_num_selected_axes'))
    num_candidate_axes = int(config.get('SymDetection-Parameters', 'intrinsic_num_candidate_axes'))
    distortion_thresh = float(config.get('SymDetection-Parameters', 'intrinsic_distortion_thresh'))
    # triplets selection
    points = '%s/fed_tips.pnts'%work_dir
    triplets = '%s/fed_triplets.triplets'%work_dir
    commands = ['%s/intrinsic_triplets_selection'%internal_bin_path ,mesh, tips,
            points, triplets,
            '-property', agd,
            '-my_triplets', '%s/my_triplets.txt'%work_dir,
            '-my_coarse', '%s/my_coarse.cor'%work_dir,
            '-v',
            '-thres_dist', min_dist_thresh,
            '-thres_dist_diff', dist_to_other_thresh,
            '-thres_prp_diff', feature_thresh,
            '-max_triplets', '10000']
    os.system(' '.join(commands))
    # intrinsic symmetry blended
    commands = ['%s/IntrinsicSymmetryBlended'%current_path, mesh, work_dir, '%s/SurfaceVis'%(external_bin_path), sym_settings_path, sym_vis_settings_path, points, triplets, min_num_shared]
    os.system(' '.join(commands))
    meshname = os.path.basename(mesh).split('.')[0]
    for i in range(num_candidate_axes):
        map_name = '%s/%s/IS_Blended_%sIS_Blended_DenseMapPreceiseAll_%d_%s_to_%s.dense.preceise.map'%(work_dir, meshname, meshname, i, meshname, meshname)
        if not os.path.exists(map_name):
            break
        sub_dir = '%s/%d'%(work_dir, i)
        os.system('mkdir -p %s'%(sub_dir))
        vertex_map_name = '%s/symmetry.map'%(sub_dir)
        sym_function = '%s/symmetry.val'%(sub_dir)
        binary_function = '%s/binary.val'%(sub_dir)
        symmetry_axis = '%s/symmetry.pid'%(sub_dir)
        # generate map
        commands = [ '%s/precisemap2vertexmap'%internal_bin_path, mesh, map_name, vertex_map_name]
        os.system(' '.join(commands))
        # map to sym_function
        commands = [ '%s/symmap2dist'%internal_bin_path, mesh, vertex_map_name, sym_function, '-max_dist %s'%map_to_dist_max ]
        os.system(' '.join(commands))
        # sym_function to binary
        commands = [ '%s/symdist2binary'%internal_bin_path, mesh, sym_function, binary_function, '-low_center' ]
        os.system(' '.join(commands))
        # sym_function to axis
        commands = ['%s/symdist2axis'%internal_bin_path, mesh, sym_function, symmetry_axis, '-extern_binary', binary_function, '-v']
        os.system(' '.join(commands))
        # check
        axis_file = open(symmetry_axis, 'r')
        lines = axis_file.read()
        if not lines.strip():
            print 'the axis is 0 length!'
            continue
        # refine axis
        commands = ['%s/axisrefine'%internal_bin_path, mesh, symmetry_axis, symmetry_axis, '-v']
        os.system(' '.join(commands))
    # rank symmetry axes based on 'confidence'
    confidence_info = '%s/%s/IS_Blended_%sIS_Blended_DenseMapPreceiseAll__%s_to_%s.dense.preceise.confidences.txt'%(work_dir, meshname, meshname, meshname, meshname)
    fin = open(confidence_info, 'r')
    lines = fin.readlines()
    fin.close()
    ## line 0: NumMaps XX
    num_obtained_axes = int(lines[0].split(' ')[1])
    if num_obtained_axes > num_candidate_axes:
        num_obtained_axes = num_candidate_axes
    ## line 1: Confidence AnotherNumber Confidence AnotherNumber ....
    words = lines[1].split(' ')
    ## parse line 1
    evaluation = []
    indices = []
    confidences = []
    axislen = []
    num_valid_axes = 0
    for i in range(num_obtained_axes):
        confidence = float(words[2*i])
        axis = '%s/%d/symmetry.pid'%(work_dir, i)
        sym = '%s/%d/symmetry.val'%(work_dir, i)
        indices.append(i)
        confidences.append(confidence)
        axislen.append(GetAxisLength(axis))
        if IsAxisValid(mesh, tips, axis, sym, '%s/%d'%(work_dir, i)):
            num_valid_axes = num_valid_axes + 1
            evaluation.append((confidence-distortion_thresh)/(1.0-distortion_thresh)*axislen[-1])
        else:
            evaluation.append(0)
    # we at least have one axis (first one)
    if num_valid_axes == 0 and num_obtained_axes >= 1:
        evaluation[0] = (confidence[0] - distortion_thresh)*float(axislen[0])
    # sort by evaluation
    indices.sort(lambda x,y: max_dec(evaluation[x], evaluation[y]))
    # select the top axes
    num_selected_axes = min(num_selected_axes, num_obtained_axes)
    suc_flag = False
    for i in range(num_selected_axes):
        if evaluation[indices[i]] <= 0.0:
            break
        suc_flag = True
        print '     Axis %d: (index: %d, confidence: %f, length: %d, eval: %f)'%(i, indices[i], confidences[indices[i]], axislen[indices[i]], evaluation[indices[i]])
        sub_result_dir = '%s/%d'%(result_dir, i)
        os.system('mkdir -p %s'%sub_result_dir)
        axis = '%s/%d/symmetry.pid'%(work_dir, indices[i])
        mapp = '%s/%d/symmetry.map'%(work_dir, indices[i])
        shutil.copy(axis, sub_result_dir)
        shutil.copy(mapp, sub_result_dir)
        outfile = open('%s/distortion.txt'%sub_result_dir, 'w')
        outfile.write(str(confidences[indices[i]]))
        outfile.close()
    return suc_flag
