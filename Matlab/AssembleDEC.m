% Generate a structure containing operators and surface info used in DEC
% calculations.

function out_struct = AssembleDEC( FaceArray, NodeArray )
%% Assemble surface data structures
    % Assemble an EdgeArray
    edges12 = FaceArray(:,[1,2]);
    edges23 = FaceArray(:,[2,3]);
    edges31 = FaceArray(:,[3,1]);
    % List of edges with redundancies
    edges_red = [edges12;edges23;edges31];
    % Sort edges and remove redundancies
    edges_red_sort = sort(edges_red,2);
    [EdgeArray, ~, i_EA] = unique(edges_red_sort,'rows');
    
    % Get edge orientations using the convention that edges point from
    % lesser node number to larger
    EdgeOrient = sign( edges_red(:,2) - edges_red(:,1) );
    EdgeOrient = reshape(EdgeOrient,[],3);
    
    % Get numbeer of nodes, edges, and faces
    num_nodes = size(NodeArray,1);
    num_edges = size(EdgeArray,1);
    num_faces = size(FaceArray,1);
    
    % Sparse matrix of nodes surrounding an edge
    I = repmat((1:num_edges)',2,1);
    J = reshape(EdgeArray,[],1);
    V = ones(size(I));
    EdgeNodes = sparse( I, J, V, num_edges, num_nodes );
    clear I J V
    
    % Sparse matrix of edges surrounding a face
    I = repmat((1:num_faces)',3,1);
    J = i_EA;
    V = ones(size(I));
    FaceEdges = sparse( I, J, V, num_faces, num_edges );
    clear I J V
    
    % Sparse matrix of nodes surrounding a face
    I = repmat((1:num_faces)',3,1);
    J = reshape(FaceArray,[],1);
    V = ones(size(I));
    FaceNodes = sparse( I, J, V, num_faces, num_nodes);
    clear I J V
    
    % Array of face unit normal vectors
    e12 = NodeArray( FaceArray(:,2), : ) - NodeArray( FaceArray(:,1), : );
    e23 = NodeArray( FaceArray(:,3), : ) - NodeArray( FaceArray(:,2), : );
    e31 = NodeArray( FaceArray(:,1), : ) - NodeArray( FaceArray(:,3), : ); 
    FaceNorms = cross(e12,-e31,2) ./ vecnorm(cross(e12,-e31,2),2,2);
    
%% Assemble boundary data structures
    % Get boundary edges by finding edges part of only one face
    b_edges = sum(FaceEdges',2) == 1; % logical array
    
    % Get boundary nodes by finding nodes on edges belonging to boundary 
    b_nodes = any(EdgeNodes(b_edges,:)',2); % logical array
    
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
        b_orient(i) = EdgeOrient(I,J);
        
        % Get an inward unit vector normal to the boundary
        b_tan = NodeArray(FaceArray(I,e_ijk{J}(2)),:) - NodeArray(FaceArray(I,e_ijk{J}(1)),:);
        b_tan = b_tan / norm(b_tan);
        b_norm(i,:) = cross(FaceNorms(I,:),b_tan);
        b_norm(i,:) = b_norm(i,:) / norm(b_norm(i,:));
    end

%% Assemble area/weight operators
    % Get the cotan at every face corner
    CotanArray = nan( num_faces, 3 );     
    CotanArray(:,1) = dot(e12,-e31,2) ./ abs(vecnorm(cross(e12,-e31,2),2,2));
    CotanArray(:,2) = dot(e23,-e12,2) ./ abs(vecnorm(cross(e23,-e12,2),2,2));
    CotanArray(:,3) = dot(e31,-e23,2) ./ abs(vecnorm(cross(e31,-e23,2),2,2));

    % Node areas on each face (and vice versa)
    I = repmat((1:num_faces)',3,1);
    J = reshape(FaceArray,[],1);
    V = [ dot(e31,e31,2).*CotanArray(:,2)/8 + dot(e12,e12,2).*CotanArray(:,3)/8;
          dot(e12,e12,2).*CotanArray(:,3)/8 + dot(e23,e23,2).*CotanArray(:,1)/8;
          dot(e23,e23,2).*CotanArray(:,1)/8 + dot(e31,e31,2).*CotanArray(:,2)/8 ];
    FaceNodeAreas = sparse( I, J, V, num_faces, num_nodes );
    NodeFaceAreas = FaceNodeAreas';
    clear I J V
    
    % Edge areas on each face (and vice versa)
    I = repmat((1:num_faces)',3,1);
    J = i_EA;
    V = [ dot(e12,e12,2).*CotanArray(:,3)/4;
          dot(e23,e23,2).*CotanArray(:,1)/4;
          dot(e31,e31,2).*CotanArray(:,2)/4 ];
    FaceEdgeAreas = sparse( I, J, V, num_faces, num_edges );
    EdgeFaceAreas = FaceEdgeAreas';
    clear I J V

    % Get the areas associated with each cell
    FaceAreas = sum( FaceNodeAreas, 2 );
    EdgeAreas = sum( EdgeFaceAreas, 2 );
    NodeAreas = sum( NodeFaceAreas, 2 );
    
    % Edge areas on each node
    NodeEdgeAreas = EdgeNodes' .* EdgeAreas' / 2;
    
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
    V = reshape(EdgeOrient,[],1);
    d1 = sparse( I, J, V, num_faces , num_edges );
    clear I J V

%% Assemble node star (neighbor) operators
    NodeStar = EdgeNodes' * EdgeNodes;
    NodeStar = NodeStar .* ~eye( num_nodes );
    
%% Assemble output structure
    out_struct = struct();
    out_struct.EdgeArray = EdgeArray;
    out_struct.EdgeVectors = EdgeVectors;
    out_struct.EdgeLengths = EdgeLengths;
    out_struct.EdgeDir = EdgeDir;
    out_struct.EdgeOrient = EdgeOrient;
    out_struct.FaceEdges = FaceEdges;
    out_struct.FaceNodes = FaceNodes;
    out_struct.EdgeNodes = EdgeNodes;
    out_struct.FaceNorms = FaceNorms;
    out_struct.b_nodes = b_nodes;
    out_struct.b_edges = b_edges;
    out_struct.b_faces = b_faces;
    out_struct.b_norm = b_norm;
    out_struct.b_orient = b_orient;
    out_struct.NodeAreas = NodeAreas;
    out_struct.EdgeAreas = EdgeAreas;
    out_struct.FaceAreas = FaceAreas;
    out_struct.FaceNodeAreas = FaceNodeAreas;
    out_struct.FaceEdgeAreas = FaceEdgeAreas;
    out_struct.NodeEdgeAreas = NodeEdgeAreas;
    out_struct.EdgeRatios = EdgeRatios;
    out_struct.hs0 = hs0;
    out_struct.hs1 = hs1;
    out_struct.hs2 = hs2;
    out_struct.d0 = d0;
    out_struct.d1 = d1;
    out_struct.NodeStar = NodeStar;
    
end