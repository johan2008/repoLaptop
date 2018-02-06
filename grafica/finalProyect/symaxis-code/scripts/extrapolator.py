from variable_def import *
from parse_work_directories import *
# Check arguments
if len(sys.argv) < 5 :
    arguments = [ 'python', sys.argv[0], 'work_folder', 'index-0-or-1', 'meshname0', 'meshname1']
    print ' '.join(arguments);
    sys.exit()

print '################################################################'
print '#                         Extrapolator                        #'
print '################################################################'

###### executables : symmap2tipscor, symaxis2tipscor, axisalign_symtips_to_tipscor, cor2triplets, ./BlendedInterpolation
# prepare variables
work_folder = sys.argv[1]
index = sys.argv[2]
meshnames = [sys.argv[3], sys.argv[4]]
[mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder] = parse_work_directories(work_folder)

results_align_folder = results_folder+'/axis_align'+index
results_map_folder = results_folder+'/map'+index
work_folder = map_work_folder+'/'+index

os.system('mkdir -p '+work_folder)
os.system('mkdir -p '+results_map_folder)

meshes = [results_folder+'/mesh0.off', results_folder+'/mesh1.off']

# Load configuration file
config = ConfigParser.RawConfigParser()
config.read('%s/config.cfg'%current_path)

# align tips
if os.path.exists(results_align_folder+'/symmetry0.map') and os.path.exists(results_align_folder+'/symmetry1.map'):
    print 'symmetry*.map is available, and we can find tips-self-cors by symmetry map!'
    for i in range(2):
        tips = results_align_folder+'/tips%d.pid'%i
        symmap = results_align_folder+'/symmetry%d.map'%i
        axis = results_align_folder+'/axis%d.pid'%i
        tips_self_cor = work_folder + '/selfcor%d.cor'%i
        arguments = ['%s/symmap2tipscor'%internal_bin_path, meshes[i], tips, symmap, axis, tips_self_cor, '-thres 0.15']
        os.system(' '.join(arguments))
else:
    print 'symmetry*.map is unavailable, and we find tips-self-cors without symmetry map!'
    for i in range(2):
        tips = results_align_folder+'/tips%d.pid'%i
        axis = results_align_folder+'/axis%d.pid'%i
        tips_self_cor = work_folder + '/selfcor%d.cor'%i
        arguments = ['%s/symaxis2tipscor'%internal_bin_path, meshes[i], tips, axis, tips_self_cor, '-thres 0.10']
        os.system(' '.join(arguments))

# create correspondences of tips between meshes
tips_align = results_map_folder + '/tips_align.cor'
arguments = ['%s/axisalign_symtips_to_tipscor'%internal_bin_path, meshes[0], meshes[1], 
             work_folder + '/selfcor0.cor', work_folder+'/selfcor1.cor',
             results_align_folder+'/align.cor', tips_align,
             '-remove_tips_on_axis',
             '-orient 0', '-v', '-debug', '-dummy_node_cost', '1.0']
os.system(' '.join(arguments))

# find triplets
num_triplets_samples = config.get('Extrapolation-Parameters', 'triplets_num_samples')
arguments = ['%s/cor2triplets'%internal_bin_path, meshes[0], meshes[1],
             results_align_folder+'/align.cor',
             results_map_folder+'/coarse.cor', work_folder+'/triplets.txt',
             '-v', '-naxissamples', num_triplets_samples, 
             '-tipscorr', tips_align, 
             '-ntriplets_each_tip', '-1', '-no_aaa','-no_tta', '-no_ttt']
os.system(' '.join(arguments))


arguments = ['%s/BlendedInterpolation'%current_path, 
        meshes[0], meshes[1], 
        results_map_folder+'/coarse.cor',
        work_folder+'/triplets.txt',
        work_folder,
        '%s/SurfaceVis'%external_bin_path,
        blended_settings_path,
        blended_vis_settings_path
        ]
#print ' '.join(arguments)
os.system(' '.join(arguments))
print 'Blended interpolation finish!'
print 'Copy ...'
finalmap = '%s/%s/BlendedInterpolated_DenseMap_%s_to_%s.dense.gaps.map'%(work_folder,
            'mesh0','mesh0','mesh1')
finalconf = '%s/%s/BlendedInterpolated_MapConfidence_TriangleArea.conf.vals'%(work_folder, 
             'mesh0')
shutil.copy(finalconf, results_map_folder+'/distortion.txt')
shutil.copy(finalmap, results_map_folder+'/final.map')
