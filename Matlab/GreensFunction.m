% Finds the nodal values of the Green's Function at the query node.

function [G, grad_G] = GreensFunction( FaceArray, NodeArray, Nq )
%% Quick settings
    % Boundary conditions to use
    % bc = 0: Potential is set to zero at one boundary nodes
    % bc = 1: Potential is set to zero at all boundary nodes 
    bc = 1;
    
%% Assemble necessary geometry info and operators
    num_nodes = size(NodeArray,1);
    num_faces = size(FaceArray,1);
    
    DEC = AssembleDEC(FaceArray,NodeArray);
    EdgeArray = DEC.EdgeArray;
    EdgeLengths = DEC.EdgeLengths;
    EdgeDir = DEC.EdgeDir;
    b_nodes = DEC.b_nodes;
    b_edges = DEC.b_edges;
    b_orient = DEC.b_orient;
    b_norm = DEC.b_norm;    
    hs1 = DEC.hs1;
    hs2 = DEC.hs2;
    d0 = DEC.d0;
    d1 = DEC.d1;    
    NodeAreas = DEC.NodeAreas;
    NodeFaceAreas = DEC.FaceNodeAreas';
    NodeFaceWeights = NodeFaceAreas ./ NodeAreas;
    NodeEdgeAreas = DEC.NodeEdgeAreas;
    
%% Compute the potential function (0-form) and its associated 1-form
    % Set boundary conditions
    G_bc = nan(num_nodes,1);
    if bitand(bc,1) == 1
        % Set potential on all boundary nodes to 0
        bc_nodes = b_nodes;
        G_bc(bc_nodes) = 0;
    else
        % Set potential to zero at a single node to fix constant fields
        bc_nodes = (1:num_nodes)' == find(b_nodes,1);
        %bc_nodes = false(num_nodes,1);
        G_bc(bc_nodes) = 0;
    end
    
    % Compute "stiffness matrix"
    K = d0' * hs1 * d0;
    
    % Compute right hand side vector (with BC's) representing a delta
    % function at the query node
    delta = zeros( num_nodes, 1 );
    delta(Nq) = 1;
    f = delta;
    f_bc = K(:,bc_nodes) * G_bc(bc_nodes);
    RHS = f - f_bc;

    % Solve for Green's Function
    G = nan(num_nodes,1);
    G(bc_nodes) = G_bc(bc_nodes);
    G(~bc_nodes) = K(~bc_nodes,~bc_nodes) \ RHS(~bc_nodes);
    
%% The associated gradient vector field
    % grad_G: Gradient of G
    % Compute the potential gradient on every face
    grad_G = zeros(num_faces,3);
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
        ga = [ G(FaceArray(i,2)) - G(FaceArray(i,1))
               G(FaceArray(i,3)) - G(FaceArray(i,1)) ];
        % Get the gradient in the u1,u2 basis
        ga = CB*ga;
        % Get the gradient in the global basis
        grad_G(i,:) = (ga(1)*u1 + ga(2)*u2)';
    end
    % diff_alpha vector field is the average of the gradients in each face
    grad_G = NodeFaceWeights * grad_G;
    
%% Alternative way of recovering gradient field
    % Get the 1-form dG
    dG = d0 * G;
    
    % Convert that 1-form to the vector componets along edges
    grad_G_edge = dG .* EdgeDir ./ EdgeLengths;
    
    % Interpolate the edge vector components to node vectors
    NodeEdgeWeights = NodeEdgeAreas ./ sum(NodeEdgeAreas,2);
%     grad_G = NodeEdgeWeights * grad_G_edge;
    
%% Test residual of Green's Function
    max_res_nbc = max( K(~bc_nodes,~bc_nodes)*G(~bc_nodes) - delta(~bc_nodes) );
    max_res = max( K*G - delta );
    
    disp(max_res)
    disp(max_res_nbc)
end