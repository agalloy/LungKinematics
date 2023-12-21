% A function that enhances the (co)exact component of HHD with
% the (co)gradient of the 1-forms squared magnitude.
% Inputs:
%   FaceArray = F x 3 nodal connectivity array
%   NodeArray = N x 3 array of nodal positions
%   omega = E x 1 array representing the input 1-form
%   alhpa = N x 1 array representing the 0-form potential
%   beta = F x 1 array representing the 2-form copotenital
%   DEC = A structure containing relevant geometric quantities and operators
% Outputs:
%   outStruct = A structure containing the following
%   alpha = The enhanced potential
%   diff_alpha = The enhanced exact field
%   beta = The enhanced copotential
%   codiff_beta = The enhanced coexact field


function [alpha_en,beta_en] = EnhanceHHD( FaceArray, NodeArray, X, omega, alpha, beta, DEC )
%% Assemble necessary geometry prerequisites
    EdgeArray = DEC.EdgeArray;
    EdgeLengths = DEC.EdgeLengths;
    EdgeDir = DEC.EdgeDir;
    EdgeVec = EdgeLengths .* EdgeDir;
    b_nodes = DEC.b_nodes;
    b_edges = DEC.b_edges;
    b_orient = DEC.b_orient;
    b_norm = DEC.b_norm;   
    hs0 = DEC.hs0;
    hs1 = DEC.hs1;
    hs2 = DEC.hs2;
    d0 = DEC.d0;
    d1 = DEC.d1;    
    NodeAreas = DEC.NodeAreas;
    EdgeAreas = DEC.EdgeAreas;
    FaceAreas = DEC.FaceAreas;
    NodeEdgeAreas = DEC.NodeEdgeAreas;
    FaceEdgeAreas = DEC.FaceEdgeAreas;
    NodeFaceAreas = DEC.FaceNodeAreas';
    
    num_nodes = size(NodeArray,1);
    num_faces = size(FaceArray,1);
    num_edges = size(EdgeArray,1);
    
%% Compute the square magnitude of input vector field at nodes and faces
    % Metric tensor on nodal vector fields
    % "Tall" vectors are vector fields reshaped from an N X 3 matrix to an 
    % 3N X 1 column vector
    g_tall = speye(num_nodes*3);
    % Get the square magnitude "m2" on nodes
    X_tall = reshape(X,[],1);
    m2_node_tall = X_tall .* g_tall * X_tall;
    m2_node = m2_node_tall(1:num_nodes) + m2_node_tall((1:num_nodes) + num_nodes) + m2_node_tall((1:num_nodes) + 2*num_nodes);
    
    % Get the magnitude on faces    
    %m2_face = FaceEdgeAreas * (m2_edge ./ EdgeAreas);
    
%% Enhance exact component
    % Remove the "DC component" from m2 and alpha
    const_field = ones(num_nodes,1) / sqrt(ones(num_nodes,1)' * hs0 * ones(num_nodes,1));
    alpha = alpha - (const_field' * hs0 * alpha) * const_field;
    m2_node = m2_node - (const_field' * hs0 * m2_node) * const_field;    
    
    % Get the gradient of the m2 and alpha fields
    diff_m2 = d0 * m2_node;  
    diff_alpha = d0 * alpha;
    
    % Get an inner product on scalar fields that is compatible with
    % gradient fields
    ip0 = d0' * hs1 * d0;
    
    % Get the projection of m2 onto alpha and its gradient
    m2_alpha = (m2_node' * ip0 * alpha) * alpha / (alpha' * ip0 * alpha);
    diff_m2_alpha = d0 * m2_alpha;
        
    % Enhance alpha
    k = ((omega' * hs1 * diff_m2) - (diff_alpha' * hs1 * diff_m2)) /...
        ((diff_m2' * hs1 * diff_m2) - (diff_m2_alpha' * hs1 * diff_m2));
    if isnan(k)
        alpha_en = alpha;
    else
        alpha_en = alpha + k * m2_node - k * m2_alpha;
    end
    
%% Enhance coexact component

%% Assemble output
    beta_en = beta;
    
%% Verification
edge_alpha = 0;

figure()   
subplot(2,4,1)
hold on
title('Magnitude squared of field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',m2_node,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

grad_m2 = GradientVectorField( FaceArray, NodeArray, m2_node, DEC );
options_sm = struct();
subplot(2,4,2)
title('Gradient of Magnitude Squared')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), grad_m2(:,[1,2]), options_sm )
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
hold off

subplot(2,4,3)
hold on
title('Potential field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

grad_alpha = GradientVectorField( FaceArray, NodeArray, alpha, DEC );
options_sm = struct();
subplot(2,4,4)
title('Gradient of Potential')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), grad_alpha(:,[1,2]), options_sm )
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
hold off

subplot(2,4,5)
hold on
title('Component of Magnitude Squared along Potential Field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',m2_alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

grad_m2_alpha = GradientVectorField( FaceArray, NodeArray,m2_alpha, DEC );
options_sm = struct();
subplot(2,4,6)
title('Gradient of Component of Magnitude Squared along Potential Field')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), grad_m2_alpha(:,[1,2]), options_sm )
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
hold off

subplot(2,4,7)
hold on
title('Component of Magnitude Squared Perp to Potential Field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',m2_node - m2_alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

grad_m2_alpha = GradientVectorField( FaceArray, NodeArray,m2_node - m2_alpha, DEC );
options_sm = struct();
subplot(2,4,8)
title('Gradient of Component of Magnitude Squared Perp to Potential Field')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), grad_m2_alpha(:,[1,2]), options_sm )
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
hold off

% % Distance of nodes from origin
% r = vecnorm( NodeArray, 2, 2 );
% figure()
% hold on
% plot( r, m2_node, 'o' )
% plot( r(~b_nodes), m2_node(~b_nodes), 'o' )
% plot( (0:0.001:1.4), (0:0.001:1.4).^2, '-' )
% hold off
% xlabel('Distance from origin')
% ylabel('Square magnitude of field')


end