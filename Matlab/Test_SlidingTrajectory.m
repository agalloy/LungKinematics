clear
clc

%% Test against lung results file
subject = 'PFS011';
side = 'Left';

mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3\';
disp_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TidalBreathing\';


inStruct = struct();
inStruct.feb_file = fullfile(mesh_dir,[subject,'_',side,'Lung_Lobes_Mesh_v3.feb']);
%inStruct.pos_data = fullfile(disp_dir,[subject,'_',side,'Lung_Lobes_f0_PositionData.txt']);
inStruct.disp_data = fullfile(disp_dir,[subject,'_',side,'Lung_Lobes_f0_FullDisp.txt']);
inStruct.ref_surface = 2;
inStruct.mov_surface = 1;
inStruct.plot_invert = false;
inStruct.plot_blank = 15;
inStruct.start_step = 1;
inStruct.two_step = false;
inStruct.projection = 'ClosestPointTriSurf';
inStruct.dist_tol = 0.5;

tic
outStruct_lung = SlidingTrajectoryV4(inStruct);
toc