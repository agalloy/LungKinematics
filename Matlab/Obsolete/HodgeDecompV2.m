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

function outStruct = HodgeDecompV2(FaceArray,NodeArray,X)
%% Quick settings
    % Boundary conditions to use
    % bc = 0: Potential is set to zero at a single boundary node 
    %   (Harmonic component is absorbed into exact component)
    % bc = 1: Potential is set to zero at all boundary nodes 
    %   (Harmonic component is orthogonal to exact and coexact components)
    bc = 1;
    
    % Display verification info
    verify = false;


%% Assemble necessary geometry prerequisites    
    DEC = AssembleDEC(FaceArray,NodeArray);
    EdgeArray = DEC.EdgeArray;
    EdgeLengths = DEC.EdgeLengths;
    EdgeDir = DEC.EdgeDir;
    EdgeVec = EdgeLengths .* EdgeDir;
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
    
    num_nodes = size(NodeArray,1);
    num_edges = size(EdgeArray,1);
    num_faces = size(FaceArray,1);
    
%% Assemble a linear operator that converts a nodal vector field to a 1-form
    % Map a (NX2) X 1 array to an E X 1 array
    I = repmat((1:num_edges)',6,1);
    J = [ 
          reshape(EdgeArray,[],1); 
          reshape(EdgeArray,[],1) + num_nodes; 
          reshape(EdgeArray,[],1) + 2*num_nodes  
         ];
    V = [ 
          repmat(EdgeVec(:,1),2,1); 
          repmat(EdgeVec(:,2),2,1)
          repmat(EdgeVec(:,3),2,1)
         ] / 2;
    VecToForm = sparse( I, J, V, num_edges, num_nodes*3 );
    clear I J V
    
    % Convert the vector field X to a 1-form
    X_tall = reshape(X,[],1);
    omega = VecToForm * X_tall;
    
    % If verifying Stoke's Theorem, get boundary flux
    if verify
        % Interpolate X to the midpoint of each edge 
        X_mid = (X(EdgeArray(:,2),:) + X(EdgeArray(:,1),:)) / 2;    
        % Break the vector into normal and tangent components
        X_n = dot( X_mid(b_edges,:), b_norm, 2 );    
        % Get the flux through each boundary edge
        omega_star = X_n .* EdgeLengths(b_edges);
    end
   
%% Compute the potential function (0-form) and its associated 1-form
    % Set boundary conditions
    alpha_bc = nan(num_nodes,1);
    if bitand(bc,1) == 1
        % Set potential on all boundary nodes to 0
        bc_nodes = b_nodes;
        alpha_bc(bc_nodes) = 0;
    else
        % Set potential to zero at a single point
        bc_nodes = (1:num_nodes)' == find(b_nodes,1);
        alpha_bc(bc_nodes) = 0;
    end
    
    % Compute "stiffness matrix"
    K = d0' * hs1 * d0;
    
    % Compute right hand side vector (with BC's)
    f = d0' * hs1 * omega;
    f_bc = K(~bc_nodes,bc_nodes) * alpha_bc(bc_nodes);
    RHS = f(~bc_nodes) - f_bc;
    
    % Solve system for alpha
    alpha = nan(num_nodes,1);
    alpha(bc_nodes) = alpha_bc(bc_nodes);
    alpha(~bc_nodes) = K(~bc_nodes,~bc_nodes) \ RHS;

    
    % Get the 1-form diff_alpha
    diff_alpha = d0 * alpha;    

%% Compute the copotential function (2-form) and its associated 1-form
    % Set boundary conditions
    beta_bc = nan(num_faces,1);
    bc_faces = false(num_faces,1);
    % No boundary conditions defined for copotential as of now
    % However, discretization automatically enforces 0 copotential along boundary 
    
    % Compute "stiffness matrix"
    K = d1 * hs1^(-1) * d1' * hs2;
    
    % Compute right hand side vector (with BC's)
    f = d1 * omega;
    f_bc = K(~bc_faces,bc_faces) * beta_bc(bc_faces);
    RHS = f(~bc_faces) - f_bc;
    
    % Solve system for beta
    beta = nan(num_faces,1);
    beta(bc_faces) = beta_bc(bc_faces);
    beta(~bc_faces) = K(~bc_faces,~bc_faces) \ RHS;

    % Get the 1-form codiff_beta
    codiff_beta = hs1^(-1) * d1' * hs2 * beta;

%% Reconstruct vector fields using the least squares method
    % Reconstruct the Exact field:
    diff_alpha_tall = lsqminnorm( VecToForm, diff_alpha, 'warn' );
    diff_alpha_v = reshape( diff_alpha_tall, [], 3 );
    
    % Reconstruct the Coexact field:
    codiff_beta_tall = lsqminnorm( VecToForm, codiff_beta, 'warn' );
    codiff_beta_v = reshape( codiff_beta_tall, [], 3 );
    
    % Harmonic field is the leftovers
    gamma = omega - diff_alpha - codiff_beta;
    gamma_v = X - diff_alpha_v - codiff_beta_v;
    
    % Nodal interpolation of copotential
    beta_n = NodeFaceWeights * hs2 * beta;
    beta_n(b_nodes) = 0;    

%% Assemble output structure
    outStruct.alpha = alpha;
    outStruct.beta = beta;
    outStruct.beta_n = beta_n;
    outStruct.gamma = gamma_v;
    outStruct.diff_alpha = diff_alpha_v;
    outStruct.codiff_beta = codiff_beta_v;
    outStruct.omega = X;
    
%% Display text outputs
    % Get an inner product operator on 1-forms
    ip1 = hs1;  
    disp('Magnitude of vector field:')
    omega_mag = sqrt( omega' * ip1 * omega );
    disp( omega_mag )
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
        
        % Check Stoke's Theorem for each vector field
        disp('Check Stokes Theorem on Omega: Curl')
        curl_omega = d1 * omega;
        disp( sum(omega(b_edges).*b_orient) - sum(curl_omega) )    
        disp('Check Stokes Theorem on Omega: Divergence')
        flux = sum( omega_star );
        div_omega = d0(~b_edges,~b_nodes)' * hs1(~b_edges,~b_edges) * omega(~b_edges);
        disp( (flux - sum(div_omega))/flux )

        % Number of boundary nodes, edges, and faces
        disp('Number of boundary nodes, edges, and faces:')
        disp( sum(b_nodes) )
        disp( sum(b_edges) )
        disp( sum(DEC.b_faces) )
        disp('Number of nodes, edges, and faces:')
        disp( num_nodes )
        disp( size(EdgeArray,1) )
        disp( num_faces )        
        
        % Inner product on vector fields
        ipv = DEC.hs0;
        omega_mag_v = sqrt( sum(dot(omega_v, ipv*omega_v, 2)) );
        disp('One-form % loss')
        disp( 100 - omega_mag / omega_mag_v * 100 )
    end
end