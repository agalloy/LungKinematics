%% Initialize Matlab
clear
clc
close all

%% User Parameters: Mesh
% Mesh bounds
mesh_bounds = [ -1, -1
                 1,  1 ];
% Number of boundary segments (for each dimension)
num_seg = [50,50];

% Set rng seed for consistency between runs 
% Used: 37484
rng(37484) 
             
%% User Parameters: Vector field
% Transformation center
O = [0,0,0];
% Uniform translation velocity
T = [1,0,0];
% Rigid rotation angular velocity (from transformation center)
R = 1;
% Scaling velocity (from transformation center)
S = [1,1,0];

% Impulse of divergence and/or curl
div_imp_mag = 0;
curl_imp_mag = 0;

% Cut off distance
d = 100;

%% User Parameters: Plot
% Seeds per dimension (must be smaller than mesh density)
num_seeds = [12,12];
% Line width
lw = 2;
% Edge Alpha
edge_alpha = 0;

%% Generate mesh
% Get node positions along sides
mx = linspace( mesh_bounds(1,1), mesh_bounds(2,1), num_seg(1) + 1 )';
my = linspace( mesh_bounds(1,2), mesh_bounds(2,2), num_seg(2) + 1 )';
% Get a rectangular grid of points
[MX,MY] = meshgrid(mx,my);
num_nodes = numel(MX);
NodeArray = [ reshape(MX,[],1), reshape(MY,[],1) ];
% Get a subset of these nodes for plotting later on
sx = round( linspace( 1, size(mx,1), num_seeds(1) )' );
sy = round( linspace( 1, size(my,1), num_seeds(2) )' );
[SX,SY] = meshgrid(sx,sy);
seed_nodes = reshape( size(my,1)*(SX-1) + SY, [], 1);
% Get a triangulated grid from the rectangular grid
FaceArray = delaunay(NodeArray);
NodeArray = [NodeArray, zeros(num_nodes,1) ];

% Jitter the interior node positions so that they no longer form right triangles
spacing = (mesh_bounds(2,:) - mesh_bounds(1,:)) ./ num_seg;
% Get boundary nodes
b_nodes = NodeArray(:,1) == mesh_bounds(1,1) | NodeArray(:,1) == mesh_bounds(2,1) ...
        | NodeArray(:,2) == mesh_bounds(1,2) | NodeArray(:,2) == mesh_bounds(2,2);
NodeArray(~b_nodes,1:2) = NodeArray(~b_nodes,1:2) + rand(size(NodeArray(~b_nodes,:),1),2) .* spacing/4;

% Resample more isotropically
ggopts = struct();
[FaceArray,NodeArray] = ggremesh(FaceArray,NodeArray,ggopts);
seed_nodes(seed_nodes > size(NodeArray,1)) = [];

%% Generate velocity field at each node
% Rigid-body rotation (uniform curl)
A_comp = R * cross( NodeArray - O, repmat([0,0,1],size(NodeArray,1),1), 2 );
% Expansion (uniform divergence)
S_comp = ( NodeArray - O ).* S;

% Divergence impulse
div_imp = div_imp_mag / (2*pi)...
          * ( NodeArray - O )...
          ./ vecnorm(NodeArray - O, 2, 2).^2;
      
% Curl impulse
curl_imp = curl_imp_mag / (2*pi)...
           * cross( NodeArray - O, repmat([0,0,1],size(NodeArray,1),1) )...
           ./ vecnorm(NodeArray - O, 2, 2).^2;

vel = A_comp + S_comp + div_imp + curl_imp + T;

% Cut off transformations at a certain distance from origin
NA_far = vecnorm(NodeArray - O, 2, 2) > d;
vel(NA_far,:) = 0;

%% Perform HHD of vector field
% Helmholtz-Hodge Decomposition
tic
%HHDStruct = HHD_GradientRecon( FaceArray, NodeArray, vel );
HHDStruct = NaturalHHD( FaceArray, NodeArray, vel );
toc

% Load results
alpha = HHDStruct.alpha;
beta = HHDStruct.beta;
beta_n = HHDStruct.beta_n;
omega = HHDStruct.omega;
diff_alpha = HHDStruct.diff_alpha;
codiff_beta = HHDStruct.codiff_beta;
gamma = HHDStruct.gamma;
% Function handle to get vector fields in MeshGrid form
VG = @(n) reshape(n,num_seg(1)+1,num_seg(1)+1);


%% View original vector field
options = struct();
options.seeds = seed_nodes;
%options.step_size = 0.001;
%options.num_steps = 3000;

% figure()
% title('Original Vector Field')
% hold on
% PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), omega(:,[1,2]), options )
% q = quiver(NodeArray(seed_nodes,1),NodeArray(seed_nodes,2),omega(seed_nodes,1),omega(seed_nodes,2),...
%     'Color','r','AutoscaleFactor',0.5);
% set(q,'LineWidth',2.5)
% hold off
% daspect([1 1 1])

%% Visualize Helmholtz-Hodge Decomposition
figure()
sp_array(1) = subplot(2,3,1);
title('Potential')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
colorbar()
    
sp_array(2) = subplot(2,3,2);
title('Original Vector Field')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), omega(:,[1,2]), options )

sp_array(3) = subplot(2,3,3);
title('Copotential')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',beta_n,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
colorbar()

sp_array(4) = subplot(2,3,4);
title('Exact vector field component (Irrotational)')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), diff_alpha(:,[1,2]), options )

sp_array(5) = subplot(2,3,5);
title('Harmonic vector field component (Incompressible and Irrotational)')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), gamma(:,[1,2]), options )

sp_array(6) = subplot(2,3,6);
title('Coexact vector field component (Incompressible)')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), codiff_beta(:,[1,2]), options )


%% Vector field with potential
% figure()
% hold on
% patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',alpha,...
%       'EdgeAlpha',edge_alpha)
% q = quiver(NodeArray(seed_nodes,1),NodeArray(seed_nodes,2),diff_alpha(seed_nodes,1),diff_alpha(seed_nodes,2),...
%     'Color','k','AutoscaleFactor',0.5);
% set(q,'LineWidth',4)
% daspect([1 1 1])
% colorbar()
%% Exact plus coexact
% figure()
% title('Exact vector field component plus Coexact vector field component')
% PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), diff_alpha(:,[1,2]) + codiff_beta(:,[1,2]), options )

%% Exact component from boundary
% tic
% close all
% A = GenerateHarmonicField( FaceArray, NodeArray, vel );
% figure()
% subplot(1,2,1)
% title('Original vector field.')
% PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), vel(:,[1,2]), options )
% subplot(1,2,2)
% title('Harmonic field generated from boundary.')
% PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), A.gamma_v(:,[1,2]), options )
% toc

%% Green's Function at a node
% xo = [0,0,0];
% Nq = dsearchn( NodeArray, xo );
% [G, grad_G] = GreensFunction( FaceArray, NodeArray, Nq );
% % Plot the Green's function with boundary conditions
% figure()   
% sp_array(1) = subplot(1,2,1);
% hold on
% title('Green''s Function with BC''s')
% p = patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',G,...
%           'EdgeAlpha',edge_alpha);
% plot3( NodeArray(Nq,1), NodeArray(Nq,2), NodeArray(Nq,3), 'r.', 'MarkerSize', 20 )
% daspect([1 1 1])
% colorbar()
% hold off
% % Plot the gradient of the Green's Function with BC's
% sp_array(2) = subplot(1,2,2);
% hold on
% title('Gradient of Green''s Function with BC''s')
% PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), grad_G(:,[1,2]), options )
% plot3( NodeArray(Nq,1), NodeArray(Nq,2), NodeArray(Nq,3), 'r.', 'MarkerSize', 20 )
% daspect([1 1 1])
% hold off
% 
% % Compute the freespace Green's Function
% G_free = -log( vecnorm( NodeArray(Nq,:) - NodeArray, 2, 2 ) )/2/pi;
% grad_G_free = GradientVectorField( FaceArray, NodeArray, G_free );
% % Plot the freespace Green's function
% figure()   
% sp_array(1) = subplot(1,2,1);
% hold on
% title('Freespace Green''s Function')
% p = patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',G_free,...
%           'EdgeAlpha',edge_alpha);
% plot3( NodeArray(Nq,1), NodeArray(Nq,2), NodeArray(Nq,3), 'r.', 'MarkerSize', 20 )
% daspect([1 1 1])
% colorbar()
% hold off
% %Plot the gradient of the freespace Green's Function with BC's
% sp_array(2) = subplot(1,2,2);
% hold on
% title('Gradient of freespace Green''s Function')
% % PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), grad_G_free(:,[1,2]), options )
% plot3( NodeArray(Nq,1), NodeArray(Nq,2), NodeArray(Nq,3), 'r.', 'MarkerSize', 20 )
% daspect([1 1 1])
% hold off
% 
% % Compare Green's function with BC's to without
% n_dist = vecnorm( NodeArray - NodeArray(Nq,:), 2, 2 );
% G_adj = G + (G_free(end) - G(end));
% 
% figure()
% hold on
% plot( n_dist, G_free, 'o' )
% plot( n_dist, G_adj, 'o' )
% legend( 'Freespace Green''s Function', 'Green''s Function with BC''s' )

%% Compare Reconstructions
% comp_link = CompareVectorReconstructions( FaceArray, NodeArray, vel );
