%% Initialize Matlab
clear
clc

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
S = [0,0,0];
% Divergence impulse
div_imp_mag = 0;
% Curl impulse
curl_imp_mag = 0;

% Use "enhanced" HHD?
enhance = true;

% Output verification info?
verify = true;

% Cut off distance
d = 100;

%% User Parameters: Plot
% Seeds per dimension (must be smaller than mesh density)
num_seeds = [11,11];
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
HHD_options = struct();
HHD_options.enhance = enhance;
HHD_options.verify = verify;
HHD_struct = HHD_GradientRecon( FaceArray, NodeArray, vel, HHD_options );
toc

% Load results
alpha = HHD_struct.alpha;
beta = HHD_struct.beta;
beta_n = HHD_struct.beta_n;
omega = HHD_struct.omega;
diff_alpha = HHD_struct.diff_alpha;
codiff_beta = HHD_struct.codiff_beta;
gamma = HHD_struct.gamma;
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
hlink = ViewHHDComponents( FaceArray, NodeArray(:,[1,2]), HHD_struct, options );
