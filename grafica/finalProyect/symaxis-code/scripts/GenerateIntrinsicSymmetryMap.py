import time,os,sys,ConfigParser
if (len(sys.argv) != 2):
    arguments = ['python', sys.argv[0], 'meshname']
    print ' '.join(arguments)
    sys.exit()
# Parameters
#    MinNumShared: this integer indicates how many generators need 
#     MinDistThresh: minimum distance between corresponding tips
#    DistToOtherThresh: threshold of corresponding geodessic distance ratio
#    FeatureThresh: threshold of agd value ratio on corresponding tips

meshname = sys.argv[1]
print '========================== Intrinsic Method : %s ==================='%(meshname)
config = ConfigParser.RawConfigParser()
configure_name = '.config.cfg'
config.read(configure_name)

models_dir = config.get('Directories', 'ModelsDir')
mesh = '%s/%s.off'%(models_dir, meshname)
work_dir = 'WorkDir/%s'%(meshname)
if os.path.exists(work_dir):
    os.system('rm -f -r %s'%(work_dir))
commands = ['mkdir', '-p', work_dir ]
os.system(' '.join(commands))

min_num_shared = config.get('Intrinsic', 'MinNumShared')
min_dist_thresh = config.get('Intrinsic', 'MinDistThresh')
dist_to_other_thresh = config.get('Intrinsic', 'DistToOtherThresh')
feature_thresh = config.get('Intrinsic', 'FeatureThresh')

#tips_dir = config.get('Directories', 'TipsFolder')
tips_dir = config.get('Directories', 'TipsFolder')
agds_dir = config.get('Directories', 'AGDsFolder')

num_axes = int(config.get('Intrinsic', 'NumAxes'))
view_flag = config.get('General', 'ViewFlag')
tips = '%s/%s.pid'%(tips_dir, meshname)
agds = '%s/%s.val'%(agds_dir, meshname)
triplet_selection_start = time.time()
# Triplets selection
commands = [ 'intrinsic_triplets_selection',mesh, tips,
             '%s/fed_tips.pnts'%(work_dir), '%s/fed_triplets.triplets'%(work_dir),
             '-property', agds,
             '-my_triplets', '%s/my_triplets.txt'%(work_dir), '-my_coarse', '%s/my_coarse.cor'%(work_dir),
             '-v', 
             '-thres_dist', min_dist_thresh, 
             '-thres_dist_diff', dist_to_other_thresh, 
             '-thres_prp_diff', feature_thresh, '-max_triplets', '10000' ]
os.system(' '.join(commands))
elapse = time.time() - triplet_selection_start;
print '***********Triplets selection time: %f seconds.'%(elapse)
# Triplets view
if view_flag == 'True':
    commands = [ 'TripletsView', mesh, mesh, '%s/my_coarse.cor'%(work_dir), '%s/my_triplets.txt'%(work_dir)]
os.system(' '.join(commands))
# Vova's script
intrinsic_symmetry_start = time.time()
commands = [ './IntrinsicSymmetryBlended', mesh, '%s/fed_tips.pnts'%(work_dir), '%s/fed_triplets.triplets'%(work_dir), min_num_shared ]
os.system(' '.join(commands))
elapse = time.time() - intrinsic_symmetry_start;
print '************IntrinsicSymmetryBlended time: %f seconds.'%(elapse)
#print ' '.join(commands)
num_candidates = int(config.get('Intrinsic', 'CandidateAxes'))
for i in range(num_candidates):
    map_name = '%s/IS_Blended_%sIS_Blended_DenseMapPreceiseAll_%d_%s_to_%s.dense.preceise.map'%(work_dir, meshname, i, meshname, meshname)
    if not os.path.exists(map_name):
        break
    sub_dir = '%s/%d'%(work_dir, i)    
    os.system('mkdir -p %s'%(sub_dir))
    vertex_map_name = '%s/symmetry.map'%(sub_dir)
    sym_function = '%s/symmetry.val'%(sub_dir)
    binary_function = '%s/binary.val'%(sub_dir)
    symmetry_axis = '%s/symmetry.pid'%(sub_dir)
    start = time.time()
    commands = [ 'precisemap2vertexmap', mesh, map_name, vertex_map_name]
    os.system(' '.join(commands))
    elapse = time.time() - start
    print '********** PreciseMap2VertexMap %d timing: %f seconds'%(i, elapse)
    # Map to sym_function
    start = time.time()
    commands = [ 'symmap2dist', mesh, vertex_map_name, sym_function, '-max_dist 0.1' ]
    os.system(' '.join(commands))
    elapse = time.time() - start
    print '********** map2sym %d timing : %f seconds'%(i, elapse)
    if view_flag == 'True':
        os.system('prpview %s %s'%(mesh, sym_function))
    # sym_function to binary_function
    start = time.time()
    commands = [ 'symdist2binary', mesh, sym_function, binary_function, '-low_center' ]
    os.system(' '.join(commands))
    if view_flag == 'True':
        os.system('prpview %s %s'%(mesh, binary_function))
    # sym_function to axis
    commands = ['symdist2axis', mesh, sym_function, symmetry_axis, '-extern_binary', binary_function, '-v']
    os.system(' '.join(commands))
    # Check whether the length of axis is 0
    axis_file = open(symmetry_axis, 'r')
    lines = axis_file.read()
    if not lines.strip():
        print 'the axis is 0 length!'
        continue
    # Refine axis
    commands = ['axisrefine', mesh, symmetry_axis, symmetry_axis, '-v']
    os.system(' '.join(commands))
    if view_flag == 'True':
        os.system('ptsview %s %s'%(mesh, symmetry_axis))
    elapse = time.time() - start
    print '************ Extracting axis %d time: %f seconds'%(i, elapse)
    # Clearing
    os.system('rm -f core.*')
