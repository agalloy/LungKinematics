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

function outStruct = HHD_GradientRecon(FaceArray,NodeArray,X,options)
%% Parse options structure or apply defaults
    % Parse optional options input
    if ~exist( 'options', 'var' )
        options = struct();
    end
    % Boundary conditions to use
        % bc = 0: Potential is set to zero at a single boundary node 
        %   (Harmonic component is absorbed into exact component)
        % bc = 1: Potential is set to zero at all boundary nodes 
        %   (Harmonic component is orthogonal to exact and coexact components)
    if isfield(options,'bc')
        bc = options.bc;
    else
        bc = 1;
    end     
    % Display verification info
    if isfield(options,'verify')
        verify = options.verify;
    else
        verify = false;
    end
    % Enhance exact and coexact fields
    if isfield(options,'enhance')
        enhance = options.enhance;
    else
        enhance = false;
    end


%% Assemble necessary geometry prerequisites
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
    
%% Convert input vector field into a discrete 1-form, omega
    % Interpolate X to the midpoint of each edge 
    X_mid = (X(EdgeArray(:,2),:) + X(EdgeArray(:,1),:)) / 2;
    
    % Break the vector into normal and tangent components
    X_t = dot( X_mid, EdgeDir, 2 );
    X_n = dot( X_mid(b_edges,:), b_norm, 2 );
    
    % Get the component of X along the edge to convert to a discrete 1-form
    omega = X_t .* EdgeLengths;
    
    % Get the flux through each boundary edge
    omega_star = X_n .* EdgeLengths(b_edges);
    
%% Perform HHD on the one-form
    HHD_options = options;
    HHD_options.DEC = DEC;
    [alpha,beta] = OneFormHHD(FaceArray,NodeArray,omega,HHD_options);
    
%% Enhance exact and coexact components
    % Enhance exact and coexact components based on the 
    % gradients of the 1-forms square magnitude field.
    if enhance
        [ alpha, diff_alpha, beta, codiff_beta ] = EnhanceHHD( FaceArray, NodeArray, X, omega, alpha, beta, DEC );
    else
        % Get the 1-forms associated with alpha and beta
        diff_alpha = d0 * alpha;  
        codiff_beta = hs1^(-1) * d1' * hs2 * beta;
    end

%% Reconstruct vector fields from (co)potentials 
    % omega: Simply use the original vector field
    omega_v = X;
    
    % grad_G: Gradient of G
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
    diff_alpha_v = NodeFaceWeights * grad_alpha;
    
    % codiff_beta: Perpendicular to gradient of beta
    % Nodal average of beta
    beta_n = NodeFaceWeights * hs2 * beta;
    %beta_n(b_nodes) = 0;
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
    codiff_beta_v = NodeFaceWeights * -cograd_beta;
    
    % gamma: Subtract the other two components from the original
    gamma = omega - diff_alpha - codiff_beta;
    gamma_v = omega_v - diff_alpha_v - codiff_beta_v;

%% Assemble output structure
    outStruct.alpha = alpha;
    outStruct.beta = beta;
    outStruct.beta_n = beta_n;
    outStruct.gamma = gamma_v;
    outStruct.gamma1F = gamma;
    outStruct.diff_alpha = diff_alpha_v;
    outStruct.diff_alpha1F = diff_alpha;
    outStruct.codiff_beta = codiff_beta_v;
    outStruct.codiff_beta1F = codiff_beta;
    outStruct.omega = omega_v;
    outStruct.omega1F = omega;
    
%% Display verification info
    % Assemble inner product operator on 1-forms
    ip1 = hs1;  
    omega_mag = sqrt( omega' * ip1 * omega );
    diff_alpha_mag = sqrt( diff_alpha' * ip1 * diff_alpha );
    codiff_beta_mag = sqrt( codiff_beta' * ip1 * codiff_beta );
    gamma_mag = sqrt( gamma' * ip1 * gamma );
    
    perc_exact = diff_alpha_mag^2 / omega_mag^2 * 100;
    perc_coexact = codiff_beta_mag^2 / omega_mag^2 * 100;
    perc_harmonic = gamma_mag^2 / omega_mag^2 * 100;
    
    if verify       
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
    
%% Display text outputs
    % Assemble inner product operator on 1-forms
    disp('Magnitude of vector field:')
    fprintf( '\t%f\n', omega_mag )
    disp('Magnitude of Exact component:')
    fprintf( '\t%f (%.1f%%)\n', diff_alpha_mag, perc_exact )
    disp('Magnitude of Coexact component:')
    fprintf( '\t%f (%.1f%%)\n', codiff_beta_mag, perc_coexact )
    disp('Magnitude of Harmonic component:')
    fprintf( '\t%f (%.1f%%)\n', gamma_mag, perc_harmonic )
end