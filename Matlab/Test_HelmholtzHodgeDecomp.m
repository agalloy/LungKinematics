% Test the Helmholtz Hodge Decomposition function
%% Initialize Matlab
clear
clc
addpath('Y:\Documents\BioMOST_Research\TecPlotTools')

%% User Parameters
% Data loading parameters
subject = 'H5983';
side = 'Left';
start_step = 10;
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3\';
disp_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TLCtoFRC_PenaltyStep\';

% Plotting parameters
plot_invert = true;

% HHD options
options.enhance = true;
options.verify = true;

% Tecplot parameters
tec_dir = 'Y:\Documents\BioMOST_Research\Lung_Analysis\Tecplot\HodgeDecomposition';
tec_pattern = '${SUBJECT}_${SIDE}Lung_t${TIME}_Enhanced';
save_tec = false;

%% Get the incremental sliding displacement field
inStruct = struct();
inStruct.feb_file = fullfile(mesh_dir,[subject,'_',side,'Lung_Lobes_Mesh_v3.feb']);
inStruct.disp_data = fullfile(disp_dir,[subject,'_',side,'Lung_Lobes_f0_FullDisp.txt']);
inStruct.ref_surface = 2;
inStruct.mov_surface = 1;
inStruct.plot_blank = 10;
inStruct.start_step = start_step;
inStruct.two_step = false;
inStruct.e_centered = false;
inStruct.plot_invert = plot_invert;
inStruct.projection = 'ClosestPointTriSurf';
inStruct.dist_tol = 0.5;
inStruct.plot_hide = true;

tic
outStruct_st = SlidingTrajectoryV4(inStruct);
toc

%% Perform Helmholtz-Hodge Decomposition
% Read sliding trajectory data
NodeArray = outStruct_st.NodeArray;
FaceArray = outStruct_st.FaceArray_Mov;
ut_inc = outStruct_st.ut_inc{1};
np_mov = outStruct_st.nodes_mov;

% Get only the fissure nodes and faces
on_fissure = all( ismember(FaceArray,np_mov), 2 );
FA_fissure = FaceArray(on_fissure,:);
% Remove nodes not attached to a face
NA_fissure = NodeArray( np_mov( ismember(np_mov,FA_fissure) ), : );
ut_inc = ut_inc( ismember(np_mov,FA_fissure), : );
np_mov = np_mov( ismember(np_mov,FA_fissure) );

% Remap global node indices to fissure indices in FA_fissure
global2local = nan( max(np_mov), 1 );
global2local(np_mov) = 1:size(np_mov,1);
FA_fissure = global2local(FA_fissure);

% Helmholtz-Hodge Decomposition
tic
HHD_struct = HHD_GradientRecon( FA_fissure, NA_fissure, ut_inc, options );
toc

alpha = HHD_struct.alpha;
beta = HHD_struct.beta;
beta_n = HHD_struct.beta_n;
omega = HHD_struct.omega;
diff_alpha = HHD_struct.diff_alpha;
codiff_beta = HHD_struct.codiff_beta;
gamma = HHD_struct.gamma;

%% Figures
options = struct();
options.spacing = 15;
%options.num_steps = 100;
options.LineWidth = 1.5;
options.plot_invert = plot_invert;
tic

hlink = ViewHHDComponents( FA_fissure, NA_fissure, HHD_struct, options );

toc

%% Extra figure
figure()
hold on
title('Original vector field')
PlotTriSurfStreamline( FA_fissure, NA_fissure, omega, options )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off

%% Compare Vector Field Reconstructions
% comp_link = CompareVectorReconstructions( FA_fissure, NA_fissure, ut_inc );

%% Save to Tecplot
if save_tec
    tec_file = replace(tec_pattern,{'${SUBJECT}','${SIDE}','${TIME}'},...
                       {subject,side,num2str(start_step)});
                   
    % Assemble data to write
    tec_struct = struct();
    tec_struct.title = tec_file;
    tec_struct.dir = tec_dir;
    tec_struct.NodeArray = NA_fissure;
    tec_struct.ElementArray = FA_fissure;
    tec_struct.VarNames = ["Potential","Copotential","Original1","Original2","Original3","Exact1","Exact2","Exact3","Coexact1","Coexact2","Coexact3","Harmonic1","Harmonic2","Harmonic3"];
    tec_struct.NodeData = [alpha,beta_n,omega,diff_alpha,codiff_beta,gamma];
    
    WriteTecPlotData(tec_struct)
end