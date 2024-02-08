% Look at the results of HHD for simple lobe transformations

%% Initialize Matlab
clear
clc
addpath('Y:\Documents\BioMOST_Research\TecPlotTools')

%% User Parameters: Data and Plotting
% Data loading parameters
subject = "H6012";
side = "Left";

% Flatten the fissure along its axis of minimal variation?
flatten = true;

% Plot parameters
plot_invert = true;

% Reference and moving surfaces:
% 1 = LTC, 2 = LLL, 3 = LUL
mov_surf = 2; % Transformations applied to this surface
ref_surf = 3; % Used to identify fissure interface

% Directories and file patterns
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3\';
mesh_pattern = '${SUBJECT}_${SIDE}Lung_Lobes_Mesh_v3.mat';

%% User Parameters: Vector field
% Transformation center ([0,0,0] is the mean of the lobe's node positions)
origin = [100,100,100];
% Scale factors (from transformation center)
scale_factors = [1,1,1];
% Rotation Axis (auto-normalized)
rot_axis = [0,0,1];
rot_axis = rot_axis / norm( rot_axis );
% Rigid rotation angular velocity (from transformation center)
rot_vel = 0.01;
% Uniform translation velocity
tran_vel = [10000,10000,10000];

% Use "enhanced" HHD?
enhance = true;

%% Load subject's mesh and identify fissure, flatten fissure if desired
% Load mesh data
mesh_name = replace( mesh_pattern, ["${SUBJECT}","${SIDE}"], [subject,side] );
mesh_file = fullfile( mesh_dir, mesh_name );
load( mesh_file, "NodeArray", "ElementArray", "nID", "eID" );

% Get the reference and moving surfaces
EA_ref = ElementArray( eID == ref_surf, : ); 
FA_ref = FESurface( EA_ref );
EA_mov = ElementArray( eID == mov_surf, : );
FA_mov = FESurface( EA_mov );

% Identify the fissure surface on the moving lobe
[ FID_fissure, NID_fissure ] = ContactInterface( FA_mov, NodeArray, FA_ref, NodeArray );
FA_fissure_global = FA_mov( FID_fissure, : );

% Get a local mesh numbering for the fissure
NA_fissure = NodeArray( NID_fissure, : );
[ ~, ~, FA_fissure_local ] = unique( FA_fissure_global );
FA_fissure_local = reshape( FA_fissure_local, [], 3 );

% Flatten fissure if desired
if flatten
    NA_fissure = FlattenPoints( NA_fissure );
end


%% Get transformation vector fields at the lobar fissure
if flatten
    center = [0,0,0];
else
    center = mean( NodeArray(nID == mov_surf,:), 1 );
end
origin_adj = origin + center;

% Get scaling vector field 
X_scale = ( NA_fissure - origin_adj ).* (scale_factors - [1,1,1]);

% Get rotation vector field
X_rot = rot_vel * cross( NA_fissure - origin_adj, repmat(rot_axis,size(NA_fissure,1),1), 2 );

% Super-impose all transformations
X = X_scale + X_rot + tran_vel;

% Restrict vector to field
X_fissure = X;

%% Perform Helmoltz-Hodge Decomposition on the generated vector field
% Perform HHD with desired options
HHD_options = struct();
HHD_options.enhance = enhance;
HHDStruct = HHD_GradientRecon( FA_fissure_local, NA_fissure, X_fissure, HHD_options );
% Read results
alpha = HHDStruct.alpha;
beta_n = HHDStruct.beta_n;
diff_alpha = HHDStruct.diff_alpha;
codiff_beta = HHDStruct.codiff_beta;
gamma = HHDStruct.gamma;

%% Plots
options = struct();
options.spacing = 15;
options.LineWidth = 1.5;

if flatten
    NA_stream = NA_fissure(:,[1,2]);
    X_stream = X_fissure(:,[1,2]);
    diff_alpha_stream = diff_alpha(:,[1,2]);
    codiff_beta_stream = codiff_beta(:,[1,2]);
    gamma_stream = gamma(:,[1,2]);
else
    NA_stream = NA_fissure;
    X_stream = X_fissure;
    diff_alpha_stream = diff_alpha;
    codiff_beta_stream = codiff_beta;
    gamma_stream = gamma;
end

figure()
title('Input vector field')
hold on
quiver3( NA_fissure(:,1), NA_fissure(:,2), NA_fissure(:,3), X_fissure(:,1), X_fissure(:,2), X_fissure(:,3) )
trimesh( FA_fissure_local, NA_fissure(:,1), NA_fissure(:,2), NA_fissure(:,3), ...
         'FaceColor', [0.80,0.80,0.80], 'EdgeColor', [0,0,0], 'FaceAlpha', 1 )
daspect([1,1,1])
xlabel('X')
ylabel('Y')
zlabel('Z')
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off

figure()   
sp_array(1) = subplot(2,3,1);
hold on
title('Potential')
p = patch('Faces',FA_fissure_local,'Vertices',NA_fissure,'FaceColor','interp','CData',alpha);
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
colorbar()
hold off

sp_array(2) = subplot(2,3,2);
hold on
title('Original vector field')
PlotTriSurfStreamline( FA_fissure_local, NA_stream, X_stream, options )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off

sp_array(3) = subplot(2,3,3);
hold on
title('Copotential')
p = patch('Faces',FA_fissure_local,'Vertices',NA_fissure,'FaceColor','interp','CData',beta_n);
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
colorbar()
hold off

sp_array(4) = subplot(2,3,4);
hold on
title('Exact vector field component (Irrotational)')
PlotTriSurfStreamline( FA_fissure_local, NA_stream, diff_alpha_stream, options )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off

sp_array(5) = subplot(2,3,5);
hold on
title('Harmonic vector field component (Incompressible and Irrotational)')
PlotTriSurfStreamline( FA_fissure_local, NA_stream, gamma_stream, options )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off

sp_array(6) = subplot(2,3,6);
hold on
title('Coexact vector field component (Incompressible)')
PlotTriSurfStreamline( FA_fissure_local, NA_stream, codiff_beta_stream, options )
daspect([1 1 1])
if plot_invert
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
end
hold off

hlink = linkprop( sp_array, {'CameraPosition','CameraUpVector'} );