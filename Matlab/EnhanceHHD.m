% A function that enhances the (co)exact component of HHD with
% the (co)gradient of the 1-forms squared magnitude.
% Inputs:
%   FaceArray = F x 3 nodal connectivity array
%   NodeArray = N x 3 array of nodal positions
%   X = N x 3 array of nodal vector field components
%   omega = E x 1 array representing the input 1-form
%   alhpa = N x 1 array representing the 0-form potential
%   beta = F x 1 array representing the 2-form copotenital
%   DEC = A structure containing relevant geometric quantities and operators
% Outputs:
%   alpha_en = The enhanced potential
%   diff_alpha_en = The enhanced exact field
%   beta_en = The enhanced copotential
%   codiff_beta_en = The enhanced coexact field


function [ alpha_en, diff_alpha_en, beta_en, codiff_beta_en ] = EnhanceHHD( FaceArray, NodeArray, X, omega, alpha, beta, DEC )
%% Assemble necessary geometry prerequisites
    b_edges = DEC.b_edges;
    b_orient = DEC.b_orient;
    hs0 = DEC.hs0;
    hs1 = DEC.hs1;
    hs2 = DEC.hs2;
    d0 = DEC.d0;
    d1 = DEC.d1;
    EdgeNodes = DEC.EdgeNodes;
    FaceEdgeAreas = DEC.FaceEdgeAreas;
    NodeAcrossEdge = DEC.NodeAcrossEdge;
    
    num_nodes = size(NodeArray,1);
    num_faces = size(FaceArray,1);
    
    HHD_opts = struct();
    HHD_opts.DEC = DEC;
    HHD_opts.bc = 1;
    
%% Compute the square magnitude of input vector field at nodes and faces
    % Get the square magnitude "m2" on nodes
    m2_node = sum( X.^2, 2 );
    
    % Interpolate vector field to circumenters and get the magnitude for each face
    FaceEdgeWeights = hs2 * FaceEdgeAreas;
    [ i, j, v ] = find( NodeAcrossEdge );
    I = i;
    J = v;
    V = FaceEdgeWeights( sub2ind(size(FaceEdgeWeights),i,j) );
    FaceNodeWeights = sparse( I, J, V, num_faces, num_nodes );
    X_face = FaceNodeWeights * X;
    m2_c_center = sum( X_face.^2, 2 );
    m2_face = hs2^(-1) * m2_c_center;
    c_center = FaceNodeWeights * NodeArray;
    
    % Interpolate vector field and magnitudes to boundary edge midpoints
    X_bedge = 1/2 * EdgeNodes(b_edges,:) * X;
    m2_bedge = sum( X_bedge.^2, 2 );
    
%% Enhance exact component 
    % Get the gradient of the m2 and alpha fields
    diff_m2 = d0 * m2_node;  
    diff_alpha = d0 * alpha;
    
    % Get the exact component of m2 with BC's
    [m2_alpha,~] = OneFormHHD( FaceArray, NodeArray, diff_m2, HHD_opts );
    diff_m2_alpha = d0 * m2_alpha;
    
    % Remove the "DC component" from potentials before adding them together
    const_field = ones(num_nodes,1) / sqrt(ones(num_nodes,1)' * hs0 * ones(num_nodes,1));
    alpha_adj = alpha - (const_field' * hs0 * alpha) * const_field;
    m2_node_adj = m2_node - (const_field' * hs0 * m2_node) * const_field;   
    m2_alpha_adj = m2_alpha - (const_field' * hs0 * m2_alpha) * const_field;
    
    % Enhance alpha
    k1 = ((omega' * hs1 * diff_m2) - (diff_alpha' * hs1 * diff_m2)) /...
        ((diff_m2' * hs1 * diff_m2) - (diff_m2_alpha' * hs1 * diff_m2));
    if isnan(k1) || isinf(k1)
        alpha_en = alpha;
        diff_alpha_en = diff_alpha;
    else
        alpha_en = alpha_adj + k1 * m2_node_adj - k1 * m2_alpha_adj;
        diff_alpha_en = d0 * alpha_en;
    end
    
%% Enhance coexact component    
    % Get the codradient of the copotential
    codiff_beta = (hs1^(-1) * d1' * hs2) * beta;
    % Get the cogradient of m2 (accounting for boundary values of m2)
    codiff_m2 = (hs1^(-1) * d1' * hs2) * m2_face;
    codiff_m2(b_edges) = codiff_m2(b_edges) - hs1(b_edges,b_edges)^(-1) * (b_orient .* m2_bedge);
    
    % Get the projection of m2 onto beta and its cogradient
    [~,m2_beta] = OneFormHHD( FaceArray, NodeArray, codiff_m2, HHD_opts );
    codiff_m2_beta = (hs1^(-1) * d1' * hs2) * m2_beta;
    
    % Remove the "DC component" from copotentials before adding them together
    const_field2 = hs2^(-1) * ones(num_faces,1);
    const_field2 = const_field2 / sqrt(const_field2' * hs2 * const_field2);
    beta_adj = beta - (const_field2' * hs2 * beta) * const_field2;
    m2_face_adj = m2_face - (const_field2' * hs2 * m2_face) * const_field2;
    m2_beta_adj = m2_beta - (const_field2' * hs2 * m2_beta) * const_field2;
    
    % Enhance beta
    k2 = ((omega' * hs1 * codiff_m2) - (codiff_beta' * hs1 * codiff_m2)) /...
        ((codiff_m2' * hs1 * codiff_m2) - (codiff_m2_beta' * hs1 * codiff_m2));
    if isnan(k2) || isinf(k2)
        beta_en = beta;
        codiff_beta_en = codiff_beta;
    else
        beta_en = beta_adj + k2 * m2_face_adj - k2 * m2_beta_adj;
        codiff_beta_en = codiff_beta + k2 * codiff_m2 - k2 * codiff_m2_beta;
    end
    
    % Check that circulation and curl are preserved with enhancement
    circ_omega = sum( b_orient.*omega(b_edges) );
    curl_omega = sum( d1 * omega );
    circ_codiff_beta = sum( b_orient.*codiff_beta(b_edges) );
    curl_codiff_beta = sum( d1 * codiff_beta );
    circ_codiff_beta_en = sum( b_orient.*codiff_beta_en(b_edges) );
    curl_codiff_beta_en = sum( d1 * codiff_beta_en );
    circ_codiff_m2 = sum( b_orient.*codiff_m2(b_edges) );
    curl_codiff_m2 = sum( d1 * codiff_m2 );
    perp_field = codiff_m2 - codiff_m2_beta;
    circ_perp_field = sum( b_orient.*perp_field(b_edges) );
    curl_perp_field = sum( d1*perp_field );
    disp( 'Curl of intermediate vector fields (first three should be equal, last should be zero):')
    disp( [curl_omega, curl_codiff_beta, curl_codiff_beta_en, curl_codiff_m2, curl_perp_field])
    disp( 'Circulation around boundary of vector fields (should be euqal to above):')
    disp( [circ_omega, circ_codiff_beta, circ_codiff_beta_en, circ_codiff_m2, circ_perp_field])
        
%% Verification
disp('Enhancement Weights (exact then coexact):')
disp([k1,k2])

edge_alpha = 0;

figure()   
subplot(2,2,1)
hold on
title('Magnitude squared of field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',m2_node,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

subplot(2,2,2)
hold on
title('Potential field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

subplot(2,2,3)
hold on
title('Component of Magnitude Squared along Potential Field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',m2_alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

subplot(2,2,4)
hold on
title('Component of Magnitude Squared Perp to Potential Field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',m2_node - m2_alpha,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off


%% Visualize Coexact enhancement
figure()   
subplot(2,2,1)
hold on
title('Magnitude squared of field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',m2_c_center,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

% subplot(2,4,2)
% hold on
% title('Circumcenter check:')
% trisurf( FaceArray, NodeArray(:,1), NodeArray(:,2), NodeArray(:,3) )
% plot3( c_center(:,1), c_center(:,2), c_center(:,3), 'r.' )
% daspect([1 1 1])
% set(gca, 'Zdir', 'reverse')
% set(gca, 'Ydir', 'reverse')

subplot(2,2,2)
hold on
title('Copotential field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',hs2*beta,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

subplot(2,2,3)
hold on
title('Component of Magnitude Squared along Copotential Field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',hs2*m2_beta,...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

subplot(2,2,4)
hold on
title('Component of Magnitude Squared Perp to Copotential Field')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',hs2*(m2_face - m2_beta),...
      'EdgeAlpha',edge_alpha);
daspect([1 1 1])
set(gca, 'Zdir', 'reverse')
set(gca, 'Ydir', 'reverse')
colorbar()
hold off

end