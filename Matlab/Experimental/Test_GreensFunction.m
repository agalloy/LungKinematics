% Visualize Green's Functions on lobar fissure surfaces.
%% Initialize Matlab
clear
clc

%% User Parameters
% Data loading parameters
subject = 'H5972';
side = 'Left';
start_step = 5;
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3\';
disp_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TLCtoFRC_PenaltyStep\';

% Plotting parameters
plot_invert = true;

% Query node
Nq = 155;

%% Get lobar fissure data
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

% Read mesh data
NodeArray = outStruct_st.NodeArray;
FaceArray = outStruct_st.FaceArray_Mov;
ut_inc = outStruct_st.ut_inc{1};
np_mov = outStruct_st.nodes_mov;

% Get only the fissure nodes and faces
on_fissure = all( ismember(FaceArray,np_mov), 2 );
FA_fissure = FaceArray(on_fissure,:);
% Remove nodes not attached to a face
NA_fissure = NodeArray( np_mov( ismember(np_mov,FA_fissure) ), : );
np_mov = np_mov( ismember(np_mov,FA_fissure) );

% Remap global node indices to fissure indices in FA_fissure
global2local = nan( max(np_mov), 1 );
global2local(np_mov) = 1:size(np_mov,1);
FA_fissure = global2local(FA_fissure);

%% Visualize to select a query node
num_nodes = size(NA_fissure,1);


figure()
hold on
title('Vertex Numbers')
p = patch('Faces',FA_fissure,'Vertices',NA_fissure,'FaceColor','interp','CData',1:num_nodes);
plot3( NA_fissure(Nq,1), NA_fissure(Nq,2), NA_fissure(Nq,3), '.' )
daspect([1 1 1])
if plot_invert
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Ydir', 'reverse')
end
colorbar()
hold off

%% Find Green's Function at query node
tic
[G, grad_G] = GreensFunction( FA_fissure, NA_fissure, Nq );
toc

%% Plot Green's Function and gradients at query node
options = struct();
options.spacing = 15;
%options.num_steps = 100;
options.LineWidth = 1.5;
tic

figure()   
sp_array(1) = subplot(1,2,1);
hold on
title('Green''s Function')
p = patch('Faces',FA_fissure,'Vertices',NA_fissure,'FaceColor','interp','CData',G);
plot3( NA_fissure(Nq,1), NA_fissure(Nq,2), NA_fissure(Nq,3), 'r.', 'MarkerSize', 20 )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
colorbar()
hold off

sp_array(2) = subplot(1,2,2);
hold on
title('Gradient of Green''s Function')
PlotTriSurfStreamline( FA_fissure, NA_fissure, grad_G, options )
plot3( NA_fissure(Nq,1), NA_fissure(Nq,2), NA_fissure(Nq,3), 'r.', 'MarkerSize', 20 )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off
