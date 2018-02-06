import os,sys,ConfigParser
# Comment Distances and Maps in the first line of settingsSymmetryBlended.txt
# So that we can only generate agds and tips
# Check arguments
if len(sys.argv) !=2:
    arguments = ['python', sys.argv[0], 'meshname']
    print ' '.join(arguments)
    sys.exit()

config = ConfigParser.RawConfigParser()
config_name = '.config.cfg'
config.read(config_name)

TipsFolder=config.get('Directories', 'TipsFolder')
AGDsFolder=config.get('Directories', 'AGDsFolder')
ModelsDir = config.get('Directories', 'ModelsDir')
# make directories in case they do not exist
os.system('mkdir -p %s'%TipsFolder)
os.system('mkdir -p %s'%AGDsFolder)

MeshName=sys.argv[1]
Mesh='%s/%s.off'%(ModelsDir, MeshName)
script_dir = 'generate_tipsagd'
WorkDir='%s/WorkDir/%s'%(script_dir,MeshName)

if os.path.exists(WorkDir):
    os.system('rm -f -r %s'%WorkDir)

os.system('mkdir -p %s'%WorkDir)

# Go to the subfolder
os.chdir(script_dir)
# run the program
arguments = ['./IntrinsicSymmetryBlended', '../%s'%Mesh]
os.system(' '.join(arguments))
# Go back to the root folder
os.chdir('..')

# copy tips (into *.pid and *.val)
agd = '%s/IS_Blended_%s.feat.vals'%(WorkDir,MeshName)
tips = '%s/IS_Blended_%s.smplset.agd'%(WorkDir,MeshName)

# convert vprp to prp
arguments = ['vprp2prp', Mesh, agd, '%s/%s.val'%(AGDsFolder, MeshName)]
os.system(' '.join(arguments))
# visualize the agd
arguments = ['prpview', Mesh, '%s/%s.val'%(AGDsFolder, MeshName)]
os.system(' '.join(arguments))
# convert vpts to pts
arguments = ['vpts2pts', Mesh, tips, '%s/%s.pid'%(TipsFolder, MeshName)]
os.system(' '.join(arguments))
# visualize the tips
arguments = ['ptsview', Mesh, '%s/%s.pid'%(TipsFolder, MeshName)]
os.system(' '.join(arguments))
