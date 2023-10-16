

function outStruct = NaturalHHD(FaceArray,NodeArray,X)
%% User inputs
    verify = true;
    
%% Assemble necessary geometry prerequisites    
    DEC = AssembleDEC(FaceArray,NodeArray);
    EdgeArray = DEC.EdgeArray;
    EdgeLengths = DEC.EdgeLengths;
    EdgeDir = DEC.EdgeDir;
    EdgeVec = EdgeLengths .* EdgeDir;
    b_nodes = DEC.b_nodes;
    b_edges = DEC.b_edges;
    b_faces = DEC.b_faces;
    b_orient = DEC.b_orient;
    b_norm = DEC.b_norm;  
    hs0 = DEC.hs0;
    hs1 = DEC.hs1;
    hs2 = DEC.hs2;
    d0 = DEC.d0;
    d1 = DEC.d1;    
    EdgeNodes = DEC.EdgeNodes;
    FaceEdges = DEC.FaceEdges;
    NodeAreas = DEC.NodeAreas;
    NodeFaceAreas = DEC.FaceNodeAreas';
    NodeFaceWeights = NodeFaceAreas ./ NodeAreas;
    
    num_nodes = size(NodeArray,1);
    num_edges = size(EdgeArray,1);
    num_faces = size(FaceArray,1);
    num_bedges = sum(b_edges);
    
    % Circumcenter dual nodes of each face
    TR = triangulation( FaceArray, NodeArray(:,1), NodeArray(:,2), NodeArray(:,3));
    dual_node = circumcenter(TR);

%% Convert the nodal vector field to a discrete 1-form
    % Create a linear operator that interpolates the nodal vector
    % field onto the edge midpoints and then takes the projection along each edge.
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
    
    % Get the edge boundary fluxes
    X_mid = (X(EdgeArray(b_edges,1),:) + X(EdgeArray(b_edges,2),:))/2;
    X_flux = dot( X_mid, -b_norm, 2 ) .* EdgeLengths(b_edges);
    
%% Compute the exact component
    % Compute the divergence of the field
    div_omega = d0' * hs1 * omega;
    % Account for boundary flux on boundary nodes
    div_omega(b_nodes) = div_omega(b_nodes) - EdgeNodes(b_edges,b_nodes)' * X_flux / 2;

    % Initialize boundary conditions
    alpha_bc = nan(num_nodes,1);  
    bc_nodes = false(num_nodes,1);
    edge_flux = zeros(num_bedges,1);
    bnode_flux = zeros(num_nodes,1);
    
    % Compute the boundary edge fluxes associated with each divergence source
    for i = 1:num_nodes
        % Get the edges of triangles connecting each edge to the source node
        TS1 = NodeArray(EdgeArray(b_edges,1),:) - NodeArray(i,:);
        TS2 = NodeArray(EdgeArray(b_edges,2),:) - NodeArray(i,:);
        % Get the angle made at the source node and triangle area
        TS_ang = acos( dot(TS1,TS2,2)...
                       ./ vecnorm(TS1,2,2)...
                       ./ vecnorm(TS2,2,2)   );
        % Assign a length-integrated flux to the boundary edge
        flux_sign = sign(dot( -b_norm, (TS1 + TS2)/2, 2 ));  
        
        % Assign a length-integrated flux to the boundary edge
        if b_nodes(i)
            i_edges = EdgeNodes(b_edges,i)==1;
            edge_flux( ~i_edges ) = edge_flux(~i_edges) + flux_sign(~i_edges) .* div_omega(i) .* TS_ang(~i_edges) / (2*pi); 
            edge_flux( i_edges ) = edge_flux(i_edges);
            % Get the angle of the remaining sector
            ext_sector = 1 - sum( flux_sign(~i_edges) .* TS_ang(~i_edges) ) / (2*pi);
            bnode_flux(i) = div_omega(i) .* ext_sector;
        else
            edge_flux = edge_flux + flux_sign .* div_omega(i) .* TS_ang / (2*pi); 
        end
    end
    
    % Set the potential at an interior node to zero to fix constant shifts
    bc_nodes( find(b_nodes,1) ) = true;
    alpha_bc( bc_nodes ) = 0;
    
    % Compute "stiffness matrix"
    K = d0' * hs1 * d0;
    
    % Assign Dirichlet BC's
    f_Dbc = K(:,bc_nodes) * alpha_bc(bc_nodes);    
    
    % Assign Neumann BC's based on edge fluxes
    node_flux = EdgeNodes(b_edges,:)' * edge_flux / 2 + bnode_flux;
    f_Nbc = node_flux;
    
    % Solve for the potential
    RHS = div_omega - f_Dbc - f_Nbc;
    alpha = nan(num_nodes,1);
    alpha(bc_nodes) = alpha_bc(bc_nodes);
    alpha(~bc_nodes) = K(~bc_nodes,~bc_nodes) \ RHS(~bc_nodes);
    
    % Exterior derivative of alpha
    diff_alpha = d0*alpha;
    
    % Solver residual
    res_p = K(~bc_nodes,~bc_nodes)*alpha(~bc_nodes) - RHS(~bc_nodes);
    
%% Compute the coexact component
    % Compute the curl of the field
    curl_omega = d1 * omega;
    
    % Initialize boundary conditions
    beta_bc = nan(num_faces,1);
    bc_faces = false(num_faces,1);
    edge_circ = zeros(num_bedges,1);
    
    % Compute the boundary edge circulations associated with each curl source
    for i = 1:num_faces
        % Get the edges of triangles connecting each edge to the source node
        TS1 = NodeArray(EdgeArray(b_edges,1),:) - dual_node(i,:);
        TS2 = NodeArray(EdgeArray(b_edges,2),:) - dual_node(i,:);
        % Get the angle made at the source node and triangle area
        TS_ang = acos( dot(TS1,TS2,2)...
                       ./ vecnorm(TS1,2,2)...
                       ./ vecnorm(TS2,2,2)   );
        % Assign a length-integrated circulation to the boundary edge
        circ_sign = sign(dot( -b_norm, (TS1 + TS2)/2, 2 )) .* b_orient;  
        
        % Assign a length-integrated circulation to the boundary edge
        edge_circ = edge_circ + circ_sign .* curl_omega(i) .* TS_ang / (2*pi); 
    end
    % Assign Dirichlet BC's to a single face to fix offset
    bc_faces(1) = true;
    beta_bc(bc_faces) = 0;
    
    % Compute "stiffness matrix"
    d1_star = d1';
    d1_star(b_edges,:) = 0;
    cd2 = hs1^(-1) * d1_star * hs2;
    K = d1 * cd2;
    
    % Assign Dirichlet BC's
    f_Dbc = K(:,bc_faces) * beta_bc(bc_faces);    
    
    % Assign Neumann BC's
    f_Nbc = d1(:,b_edges) * edge_circ;
    
    % Solve for the copotential
    RHS = curl_omega - f_Dbc - f_Nbc;
    beta = nan(num_faces,1);
    beta(bc_faces) = beta_bc(bc_faces);
    beta(~bc_faces) = K(~bc_faces,~bc_faces) \ RHS(~bc_faces);
        
    % Compute codifferential of beta
    codiff_beta = cd2 * beta;
    codiff_beta(b_edges) = edge_circ;
    
    % Solver residual
    res_c = K(~bc_faces,~bc_faces)*beta(~bc_faces) - RHS(~bc_faces);

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
    beta_n = hs0^-1 * NodeFaceAreas * hs2 * beta;
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
    outStruct.diff_alpha = diff_alpha_v;
    outStruct.codiff_beta = codiff_beta_v;
    outStruct.omega = omega_v;
    
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
        disp('Max solver residual: Potential')
        disp( max(abs( res_p )) )
        disp('Max solver residual: Copotential')
        disp( max(abs( res_c )) )
        
        % Inner product between decomposed 1-forms
        disp('Inner product: exact and coexact (should be 0 if orthogonal)')
        disp( diff_alpha' * ip1 * codiff_beta )
        disp('Inner product: exact and harmonic (should be 0 if orthogonal)')
        disp( diff_alpha' * ip1 * gamma )
        disp('Inner product: coexact and harmonic (should be 0 if orthogonal)')
        disp( codiff_beta' * ip1 * gamma )
        
        % Check Stoke's Theorem for vector fields (curl)
        net_circ = sum(omega(b_edges).*b_orient);
        disp('Check Stokes Theorem on Omega: Curl')
        curl_omega = d1 * omega;
        disp( net_circ - sum(curl_omega) )   
        disp('Check Stokes Theorem on Coexact Component: Curl')
        curl_codiff_beta = d1 * codiff_beta;
        disp( net_circ - sum(curl_codiff_beta) )
        
        % Check Stoke's Theorem for vector fields (divergence)
        net_flux = sum( X_flux );
        disp('Check Stokes Theorem on Omega: Divergence')
        disp( net_flux )
        disp( sum(div_omega) )
        disp( (net_flux + sum(div_omega)) )
        disp('Check Stokes Theorem on Exact Component: Divergence')
        div_diff_alpha = d0' * hs1 * diff_alpha + node_flux;
        disp( net_flux )
        disp( sum(node_flux) )
        disp( sum(div_diff_alpha) )
        disp( (net_flux + sum(div_diff_alpha)) )

        % Number of boundary nodes, edges, and faces
%         disp('Number of boundary nodes, edges, and faces:')
%         disp( sum(b_nodes) )
%         disp( sum(b_edges) )
%         disp( sum(b_faces) )
%         disp('Number of nodes, edges, and faces:')
%         disp( num_nodes )
%         disp( size(EdgeArray,1) )
%         disp( num_faces )        
        
        % Inner product on vector fields
        ipv = DEC.hs0;
        omega_mag_v = sqrt( sum(dot(omega_v, ipv*omega_v, 2)) );
        disp('One-form % loss')
        disp( 100 - omega_mag / omega_mag_v * 100 )
    end
    
%% Verification plots
    edge_alpha = 0.1;
    if verify
    % Input vs. output divergence    
        div_residual = div_diff_alpha - div_omega;
        figure()
        title('Divergence Residual')
        patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',div_residual,...
              'EdgeAlpha',edge_alpha);
        daspect([1,1,1])
        colorbar()
    end
end