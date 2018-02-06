import os,sys,ConfigParser,time
if (len(sys.argv) != 2):
    arguments = [ 'python', sys.argv[0], 'meshname' ]
    print ' '.join(arguments)
    sys.exit()

config = ConfigParser.RawConfigParser()
config_name = '.config.cfg'
config.read(config_name)


meshname = sys.argv[1]
print '========================== Extrinsic Method : %s ==================='%(meshname)

models_dir = config.get('Directories', 'ModelsDir')
mesh = '%s/%s.off'%(models_dir, meshname)
work_dir = 'WorkDir/%s/extrinsic'%(meshname)
commands = ['mkdir', '-p', work_dir ]
os.system(' '.join(commands))

aligned_mesh = '%s/aligned.off'%(work_dir)
ex_funs = []
ex_maps = []
ex_errors = []
ex_avg_errors = []
for i in range(3):
    sub_dir = '%s/%d'%(work_dir, i)
    os.system('mkdir -p %s'%(sub_dir))
    ex_funs.append('%s/dist.val'%(sub_dir))
    ex_maps.append('%s/symmetry.map'%(sub_dir))
    ex_errors.append('%s/error.val'%(sub_dir))
    ex_avg_errors.append('%s/avg_error.txt'%(sub_dir))
commands = ['mshalign', mesh, aligned_mesh]
os.system(' '.join(commands))
print 'Aligning mesh ok!'

commands = ['pca_extrinsic_sym_fast', aligned_mesh, '-error_properties', ex_funs[0], ex_funs[1], ex_funs[2], '-maps', ex_maps[0], ex_maps[1], ex_maps[2], '-errors', ex_errors[0], ex_errors[1], ex_errors[2], '-avg_errors', ex_avg_errors[0], ex_avg_errors[1], ex_avg_errors[2]]
os.system(' '.join(commands))
print 'Alignedmesh 2 halves ok!'

view_flag = config.get('General', 'ViewFlag')
for i in range(3):
    print 'Generating axis %d'%(i)
    sub_dir = '%s/%d'%(work_dir, i)
    start = time.time()
    binary_fun = '%s/binary.val'%(sub_dir)
    commands = ['symdist2binary', mesh, ex_funs[i], binary_fun, '-low_center']
    os.system(' '.join(commands))
    sym_axis = '%s/symmetry.pid'%(sub_dir)
    commands = ['symdist2axis', mesh, ex_funs[i], sym_axis, '-extern_binary', binary_fun]
    os.system(' '.join(commands))
    commands = ['axisrefine', mesh, sym_axis, sym_axis, '-v']
    os.system(' '.join(commands))
    infile = open(ex_avg_errors[i], 'r')
    print "Error: %s"%(infile.readline())
    elapse = time.time() - start
    print '******** Extract axis timing : %f seconds'%(elapse)
    infile.close()
