def parse_work_directories(work_folder):
    mesh_work_folder = [work_folder+'/mesh0', work_folder+'/mesh1']
    symmetry_work_folder = [work_folder+'/symmetry0', work_folder+'/symmetry1']
    axisalignment_work_folder = work_folder+'/axis_aligment'
    map_work_folder = work_folder+'/map'
    evaluation_work_folder = work_folder+'/evluation'
    results_folder = work_folder+'/results'
    return [mesh_work_folder, symmetry_work_folder, axisalignment_work_folder, map_work_folder, evaluation_work_folder, results_folder]

