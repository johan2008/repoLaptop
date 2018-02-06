from variable_def import *
# executables: pca_extrinsic_sym_fast, mshalign
def GenExtrinsicSym(mesh, work_dir, result_dir, max_error):
    # mkdir
    os.system('mkdir -p '+work_dir)
    # align mesh to axes
    aligned_mesh = '%s/aligned.off'%(work_dir)
    commands = ['%s/mshalign'%internal_bin_path, mesh, aligned_mesh]
    os.system(' '.join(commands))
    # extrinsic candidate symmetries
    axes = []
    funs = []
    binary_funs = []
    maps = []
    errors = []
    avg_errors = []
    for i in range(3):
        sub_dir = '%s/%d'%(work_dir, i)
        os.system('mkdir -p %s'%(sub_dir))
        axes.append('%s/symmetry.pid'%sub_dir)
        funs.append('%s/dist.val'%sub_dir)
        binary_funs.append('%s/binary.val'%sub_dir)
        maps.append('%s/symmetry.map'%sub_dir)
        errors.append('%s/error.val'%sub_dir)
        avg_errors.append('%s/avg_error.txt'%sub_dir)
    commands = ['%s/pca_extrinsic_sym_fast'%internal_bin_path, aligned_mesh, 
                '-error_properties', funs[0], funs[1], funs[2], 
                '-maps', maps[0], maps[1], maps[2], 
                '-errors', errors[0], errors[1], errors[2], 
                '-avg_errors', avg_errors[0], avg_errors[1], avg_errors[2]]
    os.system(' '.join(commands))
    num_valid_axes = 0
    for i in range(3):
        print 'Generating axis %d'%(i)
        commands = ['%s/symdist2binary'%internal_bin_path, mesh, funs[i], binary_funs[i], '-low_center']
        os.system(' '.join(commands))
        commands = ['%s/symdist2axis'%internal_bin_path, mesh, funs[i], axes[i], '-extern_binary', binary_funs[i]]
        os.system(' '.join(commands))
        commands = ['%s/axisrefine'%internal_bin_path, mesh, axes[i], axes[i], '-v']
        os.system(' '.join(commands))
        commands = ['%s/symmaprefine'%internal_bin_path, mesh, axes[i], maps[i], maps[i], '-v']
        os.system(' '.join(commands))
        infile = open(avg_errors[i], 'r')
        avg_error = float(infile.readline())
        infile.close()
        if avg_error < max_error:
            sub_result_dir = '%s/%d'%(result_dir, num_valid_axes)
            os.system('mkdir -p '+sub_result_dir)
            num_valid_axes = num_valid_axes + 1
            shutil.copy(axes[i], sub_result_dir)
            shutil.copy(maps[i], sub_result_dir)
            shutil.copy(avg_errors[i], sub_result_dir)
    if num_valid_axes > 0:
        return True
    else:
        return False
