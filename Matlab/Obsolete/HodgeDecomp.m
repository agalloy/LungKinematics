% Performs a Helmholtz-Hodge decomposition on a vector field defined on a
% triangulated surface.
% Inputs:
%   FaceArray = F X 3 array containing nodal connectivity info
%   NodeArray = N X 3 array containing nodal coordinates
%   X = N X 3 array of vector field values sampled at node positions
% outStruct fields:
%   alpha = The potential function defined at vertices. Shows sources of 
%       compressibility in the vector field.
%   diff_alpha = The irrotational vector field associated with alpha.
%   beta = The copotential function defined on faces. Shows sources of
%       rotation in the vector field.
%   codiff_beta = The incompressible vector field associated with beta.
%   gamma = The harmonic component of the vector field. What's left over from
%       the other two components. Translational component?

function outStruct = HodgeDecomp(FaceArray,NodeArray,X)
%% Quick settings
    % Boundary conditions to use
    % bc = 0: Potential is set to zero at a single boundary node
    % bc = 1: Potential is set to zero at all boundary nodes 
    bc = 1;
    
    % Display verification info
    verify = true;

%% Extract edge arrays and vectors
    num_nodes = size(NodeArray,1);
    num_faces = size(FaceArray,1);
    
    % Get a nodal connectivity list of edges
    % Note: Interior edges appear twice in opposite orientations
    EdgeArray = [ FaceArray(:,[1,2])
                  FaceArray(:,[2,3])
                  FaceArray(:,[3,1]) ];
              
    % Remove redundant edges and identify boundary edges
    % Start by sorting to help identify redundant edges
    [EdgeArray_sort, sort_ind] = sort( EdgeArray, 2, 'ascend' );
    % Note which edges got flipped during sorting
    edge_flip = sort_ind(:,1) == 2;
    % Get only the unique edges
    [EdgeArray, ~, ic] = unique( EdgeArray_sort , 'rows' );    
    num_edges = size( EdgeArray, 1 );
    
    % Get an array of edges that appear in each face
    FaceEdges = reshape(ic,[],3);
    % Similar array indicating which edges are flipped
    edge_flip = reshape(edge_flip,[],3);
    
    % Get vector pointing along each edge from start node to end node. 
    E = NodeArray(EdgeArray(:,2),:) - NodeArray(EdgeArray(:,1),:);
    E_mag = vecnorm( E, 2, 2 );
    E_dir = E ./ E_mag;
    
    % Get the cotangent of the interior angles of each face
    int_cotan = nan(num_faces,3);
    node_area = nan(num_faces,3);
    % Permutation operator
    eijk = [1,2,3,1,2,3];
    % Calculate the cotan of the interior angle at each face vertex
    for i = 1:3
        e12 = NodeArray(FaceArray(:,eijk(i+1)),:) - NodeArray(FaceArray(:,i),:);
        e13 = NodeArray(FaceArray(:,eijk(i+2)),:) - NodeArray(FaceArray(:,i),:);
        int_cotan(:,i) = abs(dot(e12,e13,2) ./ vecnorm( cross(e12,e13), 2, 2 ));
    end
    % Calculate the area portion of each face vertex's dual cell
    for i = 1:3
        e12 = NodeArray(FaceArray(:,eijk(i+1)),:) - NodeArray(FaceArray(:,i),:);
        e12_sq = dot(e12,e12,2);
        e13 = NodeArray(FaceArray(:,eijk(i+2)),:) - NodeArray(FaceArray(:,i),:);
        e13_sq = dot(e13,e13,2);
        cotan2 = int_cotan(:,eijk(i+1));
        cotan3 = int_cotan(:,eijk(i+2));
        node_area(:,i) = (e13_sq.*cotan3 + e12_sq.*cotan2) / 8;
    end
    
%% Find the boundary edges, faces, and nodes
    % Count the instances of each edge (regardless of orientation) in the mesh
    num_inst = histcounts(FaceEdges, (0:num_edges) + 0.5 );
    % Edges that appear only once are boundary edges
    b_edges = num_inst == 1; % Boolean array
    b_edges_ind = find(b_edges); % Index list
    
    % Any face that contains at least one boundary edge is a boundary face
    b_faces = any( ismember( FaceEdges, b_edges_ind ), 2 );
    b_faces_ind = find(b_faces);
    
    % Any node on a boundary edge is a boundary node
    b_nodes_ind = unique( EdgeArray(b_edges,:) );
    b_nodes = false( size(NodeArray,1), 1 );
    b_nodes(b_nodes_ind) = true;
        
%% Convert the vector field X defined at each node to a 1-form omega integrated over each edge
    % Interpolate X to the midpoint of each edge 
    X_mid = (X(EdgeArray(:,2),:) + X(EdgeArray(:,1),:)) / 2;
    % Get the component of X along the edge to convert to a discrete 1-form
    omega = dot( X_mid, E_dir, 2 ) .* E_mag;

%% Exterior derivative operators for 0 and 1 forms
    % 0-Form exterior derivative (num_edges x num_nodes)
    % d0 is -1 at the starting node and 1 at the end node for each edge
    d0 = sparse( num_edges, num_nodes );
    n_ind = sub2ind( size(d0), (1:num_edges)', EdgeArray(:,1) );
    d0(n_ind) = -1;
    n_ind = sub2ind( size(d0), (1:num_edges)', EdgeArray(:,2) );
    d0(n_ind) = 1;
    
    % 1-Form exterior derivative (num_faces x num_edges)
    % d1 is 1 if an edge is correctly oriented with face and -1 if not
    d1 = sparse( num_faces, num_edges );
    e_ind = sub2ind( size(d1), (1:num_faces)', FaceEdges(:,1) );
    d1(e_ind) = (-1) .^ edge_flip(:,1);
    e_ind = sub2ind( size(d1), (1:num_faces)', FaceEdges(:,2) );
    d1(e_ind) = (-1) .^ edge_flip(:,2);
    e_ind = sub2ind( size(d1), (1:num_faces)', FaceEdges(:,3) );
    d1(e_ind) = (-1) .^ edge_flip(:,3);
    
%% Hodge star operators for 0, 1, and 2 forms
    % 0-form hodge star operator
    % Diagonal matrix with the area of each node's dual cell
    hs0 = sparse( num_nodes, num_nodes);
    for i = 1:num_nodes
        % Add up area contributions from each face
        hs0(i,i) = sum(node_area(FaceArray == i));
    end

    % 1-form hodge star operator
    % Diagonal matrix with the ratio of the dual edge length to the
    % original edge length for each edge
    hs1 = sparse( num_edges, num_edges );
    % The ratio of edge lengths is found with the cotan formula:
    % |E*ij|/|Eij| = ( cot(Aj) + cot(Bj) ) / 2
    % Where Aj and Bj are the interior face angles across from the edge
    % on either side.
    for i = 1:num_edges
        % Get the cotangent of angles opposite the current edge
        EdgeFaces = any( FaceEdges==i, 2 );
        oppNodes = ~ismember( FaceArray(EdgeFaces,:), EdgeArray(i,:) );
        cotans = int_cotan(EdgeFaces,:);
        cotans = cotans(oppNodes);
        hs1(i,i) = sum(cotans)/2;        
    end
    
    % 2-form hodge star operator
    % Diagonal matrix with one over the area of the face for each face
    hs2 = sparse( num_faces, num_faces );
    % Compute area for each face
    for i = 1:num_faces
        % Get the positions of the vertices for the triangle
        FaceVerts = NodeArray( FaceArray(i,:), : );
        % Get two edges of the triangle sharing the first node
        e12 = FaceVerts(2,:) - FaceVerts(1,:);
        e13 = FaceVerts(3,:) - FaceVerts(2,:);
        % Area = 1/2 ||e12 x e13||
        area = norm( cross(e12,e13) )/2;
        % Store 1/area in hs2
        hs2(i,i) = 1 / area;
    end

%% Apply boundary conditions to the operators and omega
    if bc == 1
        % Set the potential to 0 at all boundary nodes
        bc_nodes = b_nodes;
        % Leave omega alone at all boundary edges
        bc_edges = false(num_edges,1);
        % Leave the copotential alone at all boundary faces
        % Note: the copotential is always zero immediately outside the boundary
        bc_faces = false(num_faces,1);
    elseif bc == 0
        % Set the potential to 0 at one node to fix equations
        bc_nodes = (1:num_nodes)' == b_nodes_ind(1);
        % Leave omega alone at all boundary edges
        bc_edges = false(num_edges,1);
        % Leave the copotential alone at all boundary faces
        % Note: the copotential is always zero immediately outside the boundary
        bc_faces = false(num_faces,1);
    end
    
    % Apply boundary conditions to left and right hand side vectors
    omega_bc = omega;
    omega_bc(bc_edges) = 0;
    alpha = nan(num_nodes,1);
    alpha(bc_nodes) = 0;
    beta = nan(num_faces,1);
    beta(bc_faces) = 0;

    % Apply boundary conditions to operators
    d0_bc = d0( :, ~bc_nodes );
    d1_bc = d1( ~bc_faces, : );
    hs1_bc = hs1;
    hs2_bc = hs2( ~bc_faces, ~bc_faces );
    
%% Compute the potential function (0-form) and its associated 1-form
    % Evaluate right hand side of linear system
    RHS = d0_bc' * hs1_bc * omega_bc;
    % Evaluate "Stiffness Matrix"
    K = d0_bc' * hs1_bc * d0_bc;
    % Evaluate alpha at non-boundary nodes
    alpha(~bc_nodes) = K^(-1) * RHS;
    
    % Get the 1-form associated with alpha
    diff_alpha = d0*alpha;
    
%% Compute the copotential function (2-form) and its associated 1-form
    % Evaluate the right hand side of the linear system
    RHS = d1_bc * omega_bc;
    % Evaluate "Stiffness Matrix"
    K = d1_bc * hs1_bc^(-1) * d1_bc' * hs2_bc;
    % Evaluate beta at non-boundary faces
    beta(~bc_faces) = K^(-1) * RHS;
    
    % Get the 1-form associated with beta
    codiff_beta = hs1^(-1) * d1' * hs2 * beta;
    
    % Get gamma
    gamma = omega_bc - diff_alpha - codiff_beta;

%% Reconstruct vector fields from 1-forms and (co)potentials
    % Operator for averaging face values into node values
    average_faces = hs0^(-1) * sparse( reshape(FaceArray,[],1), repmat((1:num_faces)',3,1), reshape(node_area,[],1) );    

    % omega: Simply use the original vector field
    omega_v = X;
    
    % diff_alpha: Gradient of alpha
    % Compute the potential gradient on every face
    grad_alpha = zeros(num_faces,3);
    for i = 1:num_faces
        % Get the positions of the vertices for the triangle
        FaceVerts = NodeArray( FaceArray(i,:), : );
        % Get two edges of the triangle sharing the first node
        e12 = FaceVerts(2,:) - FaceVerts(1,:);
        e13 = FaceVerts(3,:) - FaceVerts(1,:);
        % Normal vector
        N = cross(e12,e13) / norm(cross(e12,e13));
        % Orthonormal basis for the facet tangent plane
        u1 = e12 / norm(e12);
        u2 = cross(N,u1);
        % Change of basis tensor
        CB = [ e12*u1', e12*u2'
               e13*u1', e13*u2' ]^(-1);
        % Get the gradient in the e12,e13 basis
        ga = [ alpha(FaceArray(i,2)) - alpha(FaceArray(i,1))
               alpha(FaceArray(i,3)) - alpha(FaceArray(i,1)) ];
        % Get the gradient in the u1,u2 basis
        ga = CB*ga;
        % Get the gradient in the global basis
        grad_alpha(i,:) = (ga(1)*u1 + ga(2)*u2)';
    end
    % diff_alpha vector field is the average of the gradients in each face
    diff_alpha_v = average_faces * grad_alpha;
    
    % codiff_beta: Perpendicular to gradient of beta
    % Nodal average of beta
    beta_n = average_faces * hs2 * beta;
    beta_n(b_nodes) = 0;
    % Compute the copotential gradient on every face and rotate 90 degrees
    cograd_beta = zeros(num_faces,3);
    for i = 1:num_faces
        % Get the positions of the vertices for the triangle
        FaceVerts = NodeArray( FaceArray(i,:), : );
        % Get two edges of the triangle sharing the first node
        e12 = FaceVerts(2,:) - FaceVerts(1,:);
        e13 = FaceVerts(3,:) - FaceVerts(1,:);
        % Normal vector
        N = cross(e12,e13) / norm(cross(e12,e13));
        % Orthonormal basis for the facet tangent plane
        u1 = e12 / norm(e12);
        u2 = cross(N,u1);
        % Change of basis tensor
        CB = [ e12*u1', e12*u2'
               e13*u1', e13*u2' ]^(-1);
        % Get the gradient in the e12,e13 basis
        gb = [ beta_n(FaceArray(i,2)) - beta_n(FaceArray(i,1))
               beta_n(FaceArray(i,3)) - beta_n(FaceArray(i,1)) ];
        % Get the gradient in the u1,u2 basis
        gb = CB*gb;
        % Get the gradient in the global basis and rotate 90 degrees
        cograd_beta(i,:) = cross( N, (gb(1)*u1 + gb(2)*u2)' );
    end
    codiff_beta_v = average_faces * -cograd_beta;
    
    % gamma: Subtract the other two components from the original
    gamma_v = omega_v - diff_alpha_v - codiff_beta_v;

%% Assemble output structure
    outStruct.alpha = alpha;
    outStruct.beta = beta;
    outStruct.beta_n = beta_n;
    outStruct.gamma = gamma_v;
    outStruct.diff_alpha = diff_alpha_v;
    outStruct.codiff_beta = codiff_beta_v;
    outStruct.omega = omega_v;
    
%% Display inner product of the three vector field components
    % Get an inner product operator on 1-forms
    ip1 = hs1;  
    disp('Magnitude of vector field:')
    disp( sqrt( omega' * ip1 * omega ) )
    disp('Exact component of vector field:')
    disp( sqrt( diff_alpha' * ip1 * diff_alpha ) )
    disp('Coexact component of vector field')
    disp( sqrt( codiff_beta' * ip1 * codiff_beta ) )
    disp('Harmonic component of vector field')
    disp( sqrt( gamma' * ip1 * gamma ) )
    
    if verify
        % Residual for each solved quantity
        alpha_res = d0' * hs1 * omega - d0' * hs1 * d0 * alpha;
        disp('Max residual: Potential')
        disp( max( alpha_res(~bc_nodes) ) )
        beta_res = d1 * omega - d1 * hs1^(-1) * d1' * hs2 * beta;
        disp('Max residual: Copotential')
        disp( max( beta_res(~bc_faces) ) )
        
        % Inner product between decomposed 1-forms
        disp('Inner product: exact and coexact (should be 0 if orthogonal)')
        disp( diff_alpha' * ip1 * codiff_beta )
        disp('Inner product: exact and harmonic (should be 0 if orthogonal)')
        disp( diff_alpha' * ip1 * gamma )
        disp('Inner product: coexact and harmonic (should be 0 if orthogonal)')
        disp( codiff_beta' * ip1 * gamma )

    % Inner product on vector fields
%     ipv = hs0;
%     disp('Magnitude of vector field:')
%     disp( sqrt( sum(dot(omega_v, ipv*omega_v, 2)) ) )
%     disp('Exact component of vector field:')
%     disp( sqrt( sum(dot(diff_alpha_v, ipv*diff_alpha_v, 2)) ) )
%     disp('Coexact component of vector field')
%     disp( sqrt( sum(dot(codiff_beta_v, ipv*codiff_beta_v, 2)) ) )
%     disp('Harmonic component of vector field')
%     disp( sqrt( sum(dot(gamma_v, ipv*gamma_v, 2)) ) )
%     disp('Inner product: exact and coexact')
%     disp( sum(dot(diff_alpha_v, ipv*codiff_beta_v, 2)) )
%     disp('Inner product: exact and harmonic')
%     disp( sum(dot(diff_alpha_v, ipv*gamma_v, 2)) )
%     disp('Inner product: coexact and harmonic')
%     disp( sum(dot(codiff_beta_v, ipv*gamma_v, 2)) )
    
%         % Perform some sanity checks  
%         disp('Check if d o d = 0')
%         disp( max(d1*d0*alpha) )  
%         disp('Check if codiff is the adjugate of curl (two numbers should be the same)')
%         disp( (d1 * omega)' * hs2 * beta )
%         disp( (omega)' * hs1 * (hs1^(-1) * d1' * hs2 * beta) )
%         disp('Curl of exact field (should be 0):')
%         disp( max(d1 * diff_alpha) )    
%         disp('Divergence of coexact field (should be 0):')
%         disp( max(d0' * hs1 * codiff_beta) )
%         disp('Curl of harmonic field (should be 0):')
%         curl_gamma = d1 * gamma;
%         disp( max(curl_gamma) )
%         disp( max(curl_gamma(~b_faces)) )
%         disp('Divergence of harmonic field (should be 0):')
%         div_gamma = d0' * hs1 * gamma;
%         disp( max(div_gamma) )
%         disp( max(div_gamma(~b_nodes)) )
    end
end