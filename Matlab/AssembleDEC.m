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
    % Sort edges so they point from lowest node number to highest
    edges_red_sort = sort(edges_red,2);
    % Get an array of sorted edges without redundancies
    [EdgeArray, ~, i_EA] = unique(edges_red_sort,'rows');
    % i_EA maps an EdgeArray index to an edge_red index
    
    % Set EdgeOrient to 1 if the sort did not flip the edge
    % Set EdgeOrient to -1 if the sort flipped the edge
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
    CotanArray(:,1) = dot(e12,-e31,2) ./ vecnorm(cross(e12,-e31,2),2,2);
    CotanArray(:,2) = dot(e23,-e12,2) ./ vecnorm(cross(e23,-e12,2),2,2);
    CotanArray(:,3) = dot(e31,-e23,2) ./ vecnorm(cross(e31,-e23,2),2,2);

    % Node areas on each face (and vice versa)
    I = repmat((1:num_faces)',3,1);
    J = reshape(FaceArray,[],1);
    V = [ dot(e31,e31,2).*CotanArray(:,2)/8 + dot(e12,e12,2).*CotanArray(:,3)/8;
          dot(e12,e12,2).*CotanArray(:,3)/8 + dot(e23,e23,2).*CotanArray(:,1)/8;
          dot(e23,e23,2).*CotanArray(:,1)/8 + dot(e31,e31,2).*CotanArray(:,2)/8 ];
    FaceNodeAreas = sparse( I, J, V, num_faces, num_nodes );
    clear I J V
    
    % Edge areas on each face (and vice versa)
    I = repmat((1:num_faces)',3,1);
    J = i_EA;
    V = [ dot(e12,e12,2).*CotanArray(:,3)/4;
          dot(e23,e23,2).*CotanArray(:,1)/4;
          dot(e31,e31,2).*CotanArray(:,2)/4 ];
    FaceEdgeAreas = sparse( I, J, V, num_faces, num_edges );
    clear I J V

    % Get the areas associated with each cell
    FaceAreas = sum( FaceNodeAreas, 2 );
    EdgeAreas = sum( FaceEdgeAreas, 1 )';
    NodeAreas = sum( FaceNodeAreas, 1 )';
    
    % Edge areas on each node
    NodeEdgeAreas = EdgeNodes' .* EdgeAreas' / 2;
    
    % Get the edge vectors, lengths, and dual ratios
    EdgeVectors = NodeArray(EdgeArray(:,2),:) - NodeArray(EdgeArray(:,1),:);
    EdgeLengths = vecnorm(EdgeVectors,2,2);
    EdgeDir = EdgeVectors ./ EdgeLengths;
    EdgeRatios = 2*EdgeAreas ./ EdgeLengths.^2;
    
%% Assemble DEC operators
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
    
    % Hodge star for 0-forms (scalar to area)
    hs0 = speye( num_nodes ) .* NodeAreas;
    
    % Diagonal Hodge star for 1-forms (edge to dual edge)
    hs1 = speye( num_edges ) .* EdgeRatios;
    
    % Galerkin Hodge star for 1-forms (inner product on 1-forms)
    % Get inner products between 1-form edge bases integrated on given face
    ip_tri_eij_eij = nan( num_faces, 3 );
        % Inner product of e12 with e12
        ip_tri_eij_eij(:,1) = CotanArray(:,1)/12 + CotanArray(:,2)/12 + CotanArray(:,3)/4;
        % Inner product of e23 with e23
        ip_tri_eij_eij(:,2) = CotanArray(:,1)/4 + CotanArray(:,2)/12 + CotanArray(:,3)/12;
        % Inner product of e31 with e31
        ip_tri_eij_eij(:,3) = CotanArray(:,1)/12 + CotanArray(:,2)/4 + CotanArray(:,3)/12;
    % Assemble diagonal components Galerkin Hodge Star
    I_diag = (1:num_edges)';
    J_diag = (1:num_edges)';
    V_diag = zeros(num_edges,1);
    for i = 1 : (3*num_faces)
        V_diag(i_EA(i)) = V_diag(i_EA(i)) + ip_tri_eij_eij(i);
    end
    % Assemble off-diagonal components
    %EdgePairs = ( EdgeNodes * EdgeNodes' ) .* ~eye( num_edges );
    EdgePairs = ( FaceEdges' * FaceEdges ) .* ~eye( num_edges );
    [ I_off, J_off ] = find( EdgePairs );
    num_edge_pairs = size( I_off, 1 );
    V_off = nan( num_edge_pairs, 1 );
    for i = 1 : num_edge_pairs
        % Get vertices i&j from edge I (aka edge ij) and k from edge J (aka edge jk) 
        v_ij = EdgeArray(I_off(i),:);
        v_jk = EdgeArray(J_off(i),:);
        v_j = v_jk( ismember(v_jk,v_ij) ); % j is the shared vertex
        v_i = v_ij( v_ij ~= v_j ); % i belongs only to edge I
        v_k = v_jk( v_jk ~= v_j ); % k belongs only to edge J
        % Edge vectors for each pair of nodes
        e_ij = NodeArray(v_j,:) - NodeArray(v_i,:);
        e_jk = NodeArray(v_k,:) - NodeArray(v_j,:);
        e_ki = NodeArray(v_i,:) - NodeArray(v_k,:);
        % Get the cotans of the angles of the triangle spanned by edges I and J
        cot_ijk = dot( -e_ij, e_jk ) ./ vecnorm(cross( -e_ij, e_jk ));
        cot_kij = dot( -e_ki, e_ij ) ./ vecnorm(cross( -e_ki, e_ij ));
        cot_jki = dot( -e_jk, e_ki ) ./ vecnorm(cross( -e_jk, e_ki ));
        % Determine sign of the inner product:
        % Negative if they both begin/end with j, else positive
        ip_sign = (-1)^( ~xor(v_ij(1)==v_j,v_jk(1)==v_j) );
        % Inner product between edge I and edge J
        V_off(i) = ip_sign * ( cot_ijk - cot_kij - cot_jki )/12;
    end
    % Put it all together
    I = [I_diag;I_off];
    J = [J_diag;J_off];
    V = [V_diag;V_off];
    hsg = sparse( I, J, V, num_edges, num_edges );
    clear I J V 
    
    % Hodge star for 2-forms (area to scalar)
    hs2 = speye( num_faces ) ./ FaceAreas;
    
%% Assemble neighbor operators
    % Get all nodes adjacent to the selected (node star)
    NodeStar = EdgeNodes' * EdgeNodes;
    NodeStar = NodeStar .* ~eye( num_nodes );
    
    % For a given face map a node to the edge acros from it
    I = repmat((1:num_faces)',3,1);
    J1 = nan(num_faces,3);
    J2 = nan(num_faces,3);
    V1 = nan(num_faces,3);
    V2 = nan(num_faces,3);
    for i = 1:num_faces
        % Get the edges and nodes for the current face and their connectivity
        f_nodes = find( FaceNodes(i,:) )';
        f_edges = find( FaceEdges(i,:) )';
        f_EdgeNodes = EdgeNodes( f_edges, f_nodes );
        % Pair edges/nodes with the nodes/edges across from them
        J1(i,:) = f_edges;
        J2(i,:) = f_nodes;
        V1(i,:) = ~f_EdgeNodes * f_nodes;
        V2(i,:) = (~f_EdgeNodes)' * f_edges;
    end
    NodeAcrossEdge = sparse( I, J1, V1, num_faces, num_edges );
    EdgeAcrossNode = sparse( I, J2, V2, num_faces, num_nodes );
        
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
    out_struct.hsg = hsg;
    out_struct.hs2 = hs2;
    out_struct.d0 = d0;
    out_struct.d1 = d1;
    out_struct.NodeStar = NodeStar;
    out_struct.NodeAcrossEdge = NodeAcrossEdge;
    out_struct.EdgeAcrossNode = EdgeAcrossNode;
end