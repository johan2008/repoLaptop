Compling the code
    Make sure the following software and libraries are installed: OpenGL, GLUT, MATLAB, Python.
    (Skip this step for Mac OS 10.6+/Ubuntu 12.04+ users) Download and compile the source code of Surface Correspondence, and put executables under path './bin/external'.
    Type 'Make'.
    28 binaries should be generated in './bin/internal' when the code successfully compiles. 

Running the code
    Go to './scripts', and run
        python ./run_all ../data/SampleModels/386.off ../data/SampleModels/392.off ../workspace/386-392
    The first two arguments specify the source meshes, and the third argument specifies the folder where the intermediate and the final results are saved.
    After the program finishes, the final map can be found in 'workspace/386-392/results/final.map' (please look at this to understand the format of dense vertex-to-vertex map).

Running individual steps in the pipeline 
    Go to './scripts/config.cfg' and set parameters 'on' or 'off' in '[Pipeline-switches]'. Note that a latter step always relies on the results of former steps. 

Feeding groud truth at different steps
    Go to './scripts/config.cfg' and there are three modes available in '[Switches]': fully-automatic ('nothing'), perfect symmetry axes ('axis') and perfect symmetry axis alignment ('axis_alignment'). Please comment/uncomment the lines to feed ground truth data to the framework at different steps.
