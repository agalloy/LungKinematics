% Temporary function for testing basis functions produced by a set of
% boundary conditions on coexact vector fields.

function temp_CoexactBasis( FaceArray, NodeArray, Fq )
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
    
%% Compute the coexact component
    % Compute the curl of the field
    curl_omega = zeros(num_faces,1);
    curl_omega(Fq) = 1;
    
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
    bc_faces(Fq) = true;
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
    
    % Compute curl and circulation BC residuals
    curl_codiff_beta = d1*codiff_beta;
    residual = curl_codiff_beta - curl_omega;
    circ_res = codiff_beta(b_edges) - edge_circ;
    
    disp('Input curl:')
    disp( sum(curl_omega) )
    disp('Output curl:')
    disp( sum(curl_codiff_beta) )
    disp('Sum of all curl residuals:')
    disp( sum(residual) )
    disp('Max absolute curl Residual:')
    disp( max(abs(residual)) )
    disp('Max absolute circulation BC residual:')
    disp( max(abs(circ_res)) )

%% Reconstruct Vector Field
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

%% Plot functions
edge_alpha = 0.1;
options = struct();

figure()   
subplot(1,2,1);
hold on
title('Copotential basis function with circulation BC''s')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',beta_n,...
    'EdgeAlpha',edge_alpha);
plot3( dual_node(Fq,1), dual_node(Fq,2), dual_node(Fq,3), 'r.', 'MarkerSize', 20 )
daspect([1 1 1])
colorbar()
hold off
% Plot the gradient of the Green's Function with BC's
subplot(1,2,2);
hold on
title('Gradient of basis function with circulation BC''s')
PlotTriSurfStreamline( FaceArray, NodeArray(:,[1,2]), codiff_beta_v(:,[1,2]), options )
plot3( dual_node(Fq,1), dual_node(Fq,2), dual_node(Fq,3), 'r.', 'MarkerSize', 20 )
daspect([1 1 1])
hold off

% Plot distance vs. energy
beta_star = hs2 * beta;
c_dist = vecnorm( dual_node - dual_node(Fq,:), 2, 2 );
figure()
plot( c_dist, beta_star, 'o' )

figure()   
subplot(1,2,1)
hold on
title('Right hand side')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',RHS,...
    'EdgeAlpha',edge_alpha);
plot3( dual_node(Fq,1), dual_node(Fq,2), dual_node(Fq,3), 'r.', 'MarkerSize', 20 )
daspect([1 1 1])
colorbar()
hold off
subplot(1,2,2)
hold on
title('Curl Residual')
patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','flat','CData',residual,...
    'EdgeAlpha',edge_alpha);
plot3( dual_node(Fq,1), dual_node(Fq,2), dual_node(Fq,3), 'r.', 'MarkerSize', 20 )
daspect([1 1 1])
colorbar()
hold off
end