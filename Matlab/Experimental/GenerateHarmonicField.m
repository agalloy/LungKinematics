% Extend/interpolate a vector field defined on the boundary of a surface
% to a harmonic vector field within the surface

function outStruct = GenerateHarmonicField( FaceArray, NodeArray, X )
%% Assemble surface data structures
    num_nodes = size(NodeArray,1);
    num_faces = size(FaceArray,1);
    
    % Sparse matrix of nodes surrounding a face
    I = repmat((1:num_faces)',3,1);
    J = reshape(FaceArray,[],1);
    V = ones(size(I));
    FaceNodes = sparse( I, J, V, num_faces, num_nodes);
    clear I J V
    
    % Sparse matrix of faces surrounding a node
    NodeFaces = FaceNodes';
    
    % Get the number of edges and assemble an EdgeArray
    edges12 = FaceArray(:,[1,2]);
    edges23 = FaceArray(:,[2,3]);
    edges31 = FaceArray(:,[3,1]);
    % List of edges with redundancies
    edges_red = [edges12;edges23;edges31];
    % Sort edges and remove redundancies
    edges_red_sort = sort(edges_red,2);
    [EdgeArray, ~, i_EA] = unique(edges_red_sort,'rows');
    num_edges = size(EdgeArray,1);
    
    % Get edge orientations using the convention that edges point from
    % lesser node number to larger
    edge_orient = sign( edges_red(:,2) - edges_red(:,1) );
    edge_orient = reshape(edge_orient,[],3);
    
    % Sparse matrix of nodes surrounding an edge
    I = repmat((1:num_edges)',2,1);
    J = reshape(EdgeArray,[],1);
    V = ones(size(I));
    EdgeNodes = sparse( I, J, V, num_edges, num_nodes );
    clear I J V
    
    % Sparse matrix of edges surrounding a node
    NodeEdges = EdgeNodes';
    
    % Sparse matrix of edges surrounding a face
    I = repmat((1:num_faces)',3,1);
    J = i_EA;
    V = ones(size(I));
    FaceEdges = sparse( I, J, V, num_faces, num_edges );
    clear I J V
    
    % Sparse matrix of faces surrounding an edge
    EdgeFaces = FaceEdges';
    
    % Array of face unit normal vectors
    e12 = NodeArray( FaceArray(:,2), : ) - NodeArray( FaceArray(:,1), : );
    e23 = NodeArray( FaceArray(:,3), : ) - NodeArray( FaceArray(:,2), : );
    e31 = NodeArray( FaceArray(:,1), : ) - NodeArray( FaceArray(:,3), : ); 
    FaceNorms = cross(e12,-e31,2) ./ vecnorm(cross(e12,-e31,2),2,2);
    
%% Assemble boundary data structures
    % Get boundary edges by finding edges part of only one face
    b_edges = sum(EdgeFaces,2) == 1; % logical array
    
    % Get boundary nodes by finding nodes on edges belonging to boundary 
    b_nodes = any(NodeEdges(:,b_edges),2); % logical array
    
    % Get boundary faces by finding faces with edges belonging to boundary
    b_faces = any(FaceEdges(:,b_edges),2); % logical array
    
    % Get the normal and orientation of boundary edges
    b_edges_ind = find(b_edges);
    b_norm = zeros( size(b_edges_ind,1), 3 );
    b_orient = zeros( size(b_edges_ind,1), 1 );
    % Side/node permutations
    e_ijk = {[1,2],[2,3],[3,1],[1,2]};
    for i = 1:size(b_edges_ind,1)
        % Current edge id
        edge_id = b_edges_ind(i);
        % Get the side of the face associated with the edge
        [I, J] = ind2sub( [num_faces,3], find(i_EA == edge_id) );
        
        % Edge orientation from edge_flip
        b_orient(i) = edge_orient(I,J);
        
        % Get an inward unit vector normal to the boundary
        b_tan = NodeArray(FaceArray(I,e_ijk{J}(2)),:) - NodeArray(FaceArray(I,e_ijk{J}(1)),:);
        b_tan = b_tan / norm(b_tan);
        b_norm(i,:) = cross(FaceNorms(I,:),b_tan);
        b_norm(i,:) = b_norm(i,:) / norm(b_norm(i,:));
    end
    
%% Assemble area/weight operators
    % Get the cotan at every face corner
    cotan_array = nan( num_faces, 3 );     
    cotan_array(:,1) = dot(e12,-e31,2) ./ abs(vecnorm(cross(e12,-e31,2),2,2));
    cotan_array(:,2) = dot(e23,-e12,2) ./ abs(vecnorm(cross(e23,-e12,2),2,2));
    cotan_array(:,3) = dot(e31,-e23,2) ./ abs(vecnorm(cross(e31,-e23,2),2,2));

    % Node weights on each face (and vice versa)
    I = repmat((1:num_faces)',3,1);
    J = reshape(FaceArray,[],1);
    V = [ dot(e31,e31,2).*cotan_array(:,2)/8 + dot(e12,e12,2).*cotan_array(:,3)/8;
          dot(e12,e12,2).*cotan_array(:,3)/8 + dot(e23,e23,2).*cotan_array(:,1)/8;
          dot(e23,e23,2).*cotan_array(:,1)/8 + dot(e31,e31,2).*cotan_array(:,2)/8 ];
    FaceNodeAreas = sparse( I, J, V, num_faces, num_nodes );
    NodeFaceAreas = FaceNodeAreas';
    clear I J V
    
    % Edge weights on each face (and vice versa)
    I = repmat((1:num_faces)',3,1);
    J = i_EA;
    V = [ dot(e12,e12,2).*cotan_array(:,3)/4;
          dot(e23,e23,2).*cotan_array(:,1)/4;
          dot(e31,e31,2).*cotan_array(:,2)/4 ];
    FaceEdgeAreas = sparse( I, J, V, num_faces, num_edges );
    EdgeFaceAreas = FaceEdgeAreas';

    % Get the areas associated with each cell
    FaceAreas = sum( FaceNodeAreas, 2 );
    EdgeAreas = sum( EdgeFaceAreas, 2 );
    NodeAreas = sum( NodeFaceAreas, 2 );
    
    % Get the edge vectors, lengths, and dual ratios
    EdgeVectors = NodeArray(EdgeArray(:,2),:) - NodeArray(EdgeArray(:,1),:);
    EdgeLengths = vecnorm(EdgeVectors,2,2);
    EdgeDir = EdgeVectors ./ EdgeLengths;
    EdgeRatios = 2*EdgeAreas ./ EdgeLengths.^2;

%% Assemble DEC operators
    % Hodge star for 0-forms (scalar to area)
    hs0 = speye( num_nodes ) .* NodeAreas;
    
    % Hodge star for 1-forms (edge to dual edge)
    hs1 = speye( num_edges ) .* EdgeRatios;
    
    % Hodge star for 2-forms (area to scalar)
    hs2 = speye( num_faces ) ./ FaceAreas;
    
    % Exterior derivative for 0-forms
    % Subtract value at start node of edge from value at end node
    % EdgeArray gives orientation for edges
    I = repmat((1:num_edges)',2,1);
    J = reshape(EdgeArray,[],1);
    V = [ -1*ones(num_edges,1); ones(num_edges,1) ];
    d0 = sparse( I, J, V, num_edges, num_nodes );
    clear I J V
    
    % Exterior derivative for 1-forms
    % Sum the oriented 1-form values around a face
    % edge_flip stores face orientation info
    I = repmat((1:num_faces)',3,1);
    J = (1:num_edges)';
    J = J(i_EA);
    V = reshape(edge_orient,[],1);
    d1 = sparse( I, J, V, num_faces , num_edges );
    clear I J V

%% Find a unique harmonic vector field from the vectors on the boundary
    % Interpolate X to the midpoint of each edge 
    X_mid = (X(EdgeArray(:,2),:) + X(EdgeArray(:,1),:)) / 2;
    % Break the vector into normal and tangent components
    X_t = dot( X_mid, EdgeDir, 2 );
    X_n = dot( X_mid(b_edges,:), b_norm, 2 );
    % Approximate the curl on each boundary edge midpoint
    X_diff_mid = (X(EdgeArray(b_edges,2),:) - X(EdgeArray(b_edges,1),:)) ./ EdgeLengths(b_edges);
    b_curl = dot( X_diff_mid, b_norm, 2 );
    
    % Get the component of X along the edge to convert to a discrete 1-form
    omega = X_t .* EdgeLengths;
    omega_star = X_n .* EdgeLengths(b_edges);
    
    % Boundary conditions on operators: Add flux terms to d1_star
    OmegaToOmegaStar = omega_star ./ (hs1(b_edges,b_edges)*omega(b_edges));
    d1_star = d0';
    d1_star(b_nodes,b_edges) = d1_star(b_nodes,b_edges)...
                               + abs(d1_star(b_nodes,b_edges)).*OmegaToOmegaStar'/2;
    
    % Boundary conditions on operators: Add boundary curl terms to d0_star
    CurlToBCurl = (d1(b_faces,b_edges)*b_curl) ./ (hs2(b_faces,b_faces)*d1(b_faces,b_edges)*omega(b_edges));
    d0_star = d1';
    d0_star(b_edges,b_faces) = d0_star(b_edges,b_faces)...
                               + abs(d0_star(b_edges,b_faces)).*CurlToBCurl;
                           
    % Laplacian Operator
    lap = (d0 * hs0^-1 * d1_star * hs1) + (hs1^-1 * d0_star * hs2 * d1);
    
    % Evaluate right hand side vector
    BCs = lap( :, b_edges ) * omega(b_edges);
    RHS = zeros( num_edges, 1 ) - BCs;
    
    % Interpolate the boundary values to make a harmonic 1-form
    gamma = nan(size(omega));
    gamma(b_edges) = omega(b_edges);
    gamma(~b_edges) = lap(~b_edges,~b_edges) \ RHS(~b_edges);
    
%% Convert the one-form back into a vector field
    % Convert gamma to a vector field defined on edge midpoints
    gamma_ve = gamma .* EdgeDir ./ EdgeLengths;
    
    % Get the area overlap of each node and edge
    NodeEdgeAreas = NodeEdges .* EdgeAreas' / 2;
    NodeEdgeWeights = NodeEdgeAreas ./ sum(NodeEdgeAreas,2);
    % Interpolate the edge vectors to node vectors
    gamma_v = NodeEdgeWeights * gamma_ve;
    
%% Assemble output structure
    outStruct = struct();
    outStruct.gamma = gamma;
    outStruct.gamma_v = gamma_v;

%% Verification and debugging
%     disp('Check if face areas agree.')
%     FaceAreas1 = sum(FaceNodeAreas,2);
%     FaceAreas2 = sum(FaceEdgeAreas,2);
%     disp( max(abs( FaceAreas2 - FaceAreas1 )) )
%     disp( [min(FaceAreas1),max(FaceAreas1)] )
%     disp( [min(FaceAreas2),max(FaceAreas2)] )
%     
%     disp('Check if node areas agree.')
%     NodeAreas1 = sum(NodeFaceAreas,2);
%     NodeAreas2 = sum(NodeEdgeAreas,2);
%     disp( max(abs( NodeAreas2 - NodeAreas1 )) )
%     disp( [min(NodeAreas1),max(NodeAreas1)] )
%     disp( [min(NodeAreas2),max(NodeAreas2)] )
    
    disp('Check that the solution satisfies equations.')
    r_lap = lap * gamma;
    disp(max( r_lap ))
    disp(max( r_lap(~b_edges) ))
    
    disp('Check that the divergence is 0')
    div = hs0^-1 * d1_star * hs1 * gamma;
    disp(max( div ))
    disp(max( div(~b_nodes) ))
    
    disp('Check that the curl is 0')
    curl = d1 * gamma;
    disp(max( curl ))
    disp(max( curl(~b_faces) ))
    
    disp('Check Stokes Theorem: Curl')
    curl_omega = d1 * omega;
    disp( sum(omega(b_edges).*b_orient) - sum(curl_omega) )
    
    disp('Check Stokes Theorem: Divergence')
    div_omega = d1_star * hs1 * omega;
    disp( sum(div_omega) - sum(omega_star) )
    
    disp('Laplacian of Original Field')
    lap_omega = lap*omega;
    disp(max( lap_omega(b_edges) ))
    disp(max( lap_omega(~b_edges) ))
    
    disp('Max Divergence of Original Field')
    disp(max( div_omega(b_nodes) ))
    disp(max( div_omega(~b_nodes) ))
    
    disp('Max Curl of Original Field')
    disp(max( curl_omega(b_faces) ))
    disp(max( curl_omega(~b_faces) ))
       
%% Create plots of divergence and curl
    figure()
    subplot(1,2,1);
    title('Divergence')
    patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',div,...
      'EdgeAlpha',0.2);
    daspect([1 1 1])
    colorbar()
    
    subplot(1,2,2);
    title('Curl')
    patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',hs2*curl,...
      'EdgeAlpha',0.2);
    daspect([1 1 1])
    colorbar()

end