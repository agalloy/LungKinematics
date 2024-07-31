% Finds the nodal values of the Green's Function at the query node.

function [G, grad_G] = GreensFunction( FaceArray, NodeArray, Nq )
%% Quick settings
    % Boundary conditions to use
    % bc = 0: Potential is set to zero at one boundary nodes
    % bc = 1: Potential is set to zero at all boundary nodes 
    % bc = 2: Potential is set to produce radially symmetric functions
    %         (bc = 2 Only works for planar domains)
    bc = 2;
    
%% Assemble necessary geometry info and operators
    DEC = AssembleDEC(FaceArray,NodeArray);
    EdgeArray = DEC.EdgeArray;
    EdgeLengths = DEC.EdgeLengths;
    EdgeDir = DEC.EdgeDir;
    b_nodes = DEC.b_nodes;
    b_edges = DEC.b_edges;
    b_orient = DEC.b_orient;
    b_norm = DEC.b_norm; 
    hs0 = DEC.hs0;
    hs1 = DEC.hs1;
    hs2 = DEC.hs2;
    d0 = DEC.d0;
    d1 = DEC.d1;    
    EdgeNodes = DEC.EdgeNodes;
    NodeAreas = DEC.NodeAreas;
    NodeFaceAreas = DEC.FaceNodeAreas';
    NodeFaceWeights = NodeFaceAreas ./ NodeAreas;
    NodeEdgeAreas = DEC.NodeEdgeAreas;
    NodeStar = DEC.NodeStar;
    
    num_nodes = size(NodeArray,1);
    num_edges = size(EdgeArray,1);
    num_faces = size(FaceArray,1);
    num_bedges = sum(b_edges);
    
%% Compute the potential function (0-form) and its associated 1-form
    % Set boundary conditions
    G_bc = nan(num_nodes,1);  
    bc_nodes = false(num_nodes,1);
    edge_flux = zeros(num_bedges,1);
    if bitand(bc,1) == 1
        % Set potential on all boundary nodes to 0
        bc_nodes(b_nodes) = true;
        G_bc(bc_nodes) = 0;
    elseif bitand(bc,2) == 2
        % Set the potential on Nq to zero to fix constant shifts
        bc_nodes(1) = true;
        G_bc(1) = 0;
        
        % Set fluxes along each boundary edge to distribute flux with radial symmetry
        % Get the edges of triangles connecting each edge to the source node
        TS1 = NodeArray(EdgeArray(b_edges,1),:) - NodeArray(Nq,:);
        TS2 = NodeArray(EdgeArray(b_edges,2),:) - NodeArray(Nq,:);
        % Get the angle made at the source node and triangle area
        TS_ang = acos( dot(TS1,TS2,2)...
                       ./ vecnorm(TS1,2,2)...
                       ./ vecnorm(TS2,2,2)   );
        % Assign a length-integrated flux to the boundary edge
        flux_sign = sign(dot( -b_norm, (TS1 + TS2)/2, 2 ));
        edge_flux = flux_sign .* TS_ang / (2*pi);          
    elseif bitand(bc,4) == 4
        % Set the potential on Nq to zero
        %bc_nodes(Nq) = true;
        %G_bc(Nq) = 0;
        % Set the gradient from Nq to its neighbors to non-zero
        star_nodes = NodeStar(:,Nq) == 1;
        bc_edges = EdgeNodes( :, Nq ) == 1;
        dG_star_bc = -NodeEdgeAreas(Nq,bc_edges)' / NodeAreas(Nq);
        dG_bc = hs1(bc_edges,bc_edges)^(-1) * dG_star_bc;
        % Apply gradient to neighboring nodes
        bc_nodes( star_nodes ) = true;
        G_bc( star_nodes ) = dG_bc + 0;
    else
        % Set potential to zero at a single node to fix constant fields
        bc_nodes( find(b_nodes,1) ) = true;
        %bc_nodes = false(num_nodes,1);
        G_bc(bc_nodes) = 0;
    end
    
    % Compute "stiffness matrix"
    K = d0' * hs1 * d0;
    
    % Compute right hand side vector representing a delta function at the query node
    delta = zeros( num_nodes, 1 );
    delta(Nq) = 1;
    f = delta;
    
    % Assign Dirichlet BC's
    f_Dbc = K(:,bc_nodes) * G_bc(bc_nodes);    
    
    % Assign Neumann BC's based on edge fluxes
    f_Nbc = EdgeNodes(b_edges,:)' * edge_flux / 2;
    
    % Solve for Green's Function
    RHS = f - f_Dbc - f_Nbc;
    G = nan(num_nodes,1);
    G(bc_nodes) = G_bc(bc_nodes);
    G(~bc_nodes) = K(~bc_nodes,~bc_nodes) \ RHS(~bc_nodes);
    
    % Compute the gradient of Green's Function
    grad_G = GradientVectorField( FaceArray, NodeArray, G, DEC );
        
%% Alternative way of recovering gradient field
%     % Get the 1-form dG
%     dG = d0 * G;
%     
%     % Convert that 1-form to the vector componets along edges
%     grad_G_edge = dG .* EdgeDir ./ EdgeLengths;
%     
%     % Interpolate the edge vector components to node vectors
%     NodeEdgeWeights = NodeEdgeAreas ./ sum(NodeEdgeAreas,2);
%     grad_G = NodeEdgeWeights * grad_G_edge;
    
%% Test residual of Green's Function
    divG = K*G;
    res = K*G - RHS;
    max_res = max( abs(res) );
    
    disp('Maximum absolute Green''s function residual:')
    disp( abs(max_res) )
    disp('Divergence at source node')
    disp( divG(Nq) )
    disp('Sum of boundary fluxes')
    disp( sum(edge_flux) )
    disp('Sum of RHS interior')
    disp( sum(RHS(~b_nodes)) )
    
%% Verification plots
%     edge_alpha = 0.1;
%     figure()   
%     hold on
%     title('Divergence Residual')
%     patch('Faces',FaceArray,'Vertices',NodeArray,'FaceColor','interp','CData',res,...
%               'EdgeAlpha',edge_alpha);
%     plot3( NodeArray(Nq,1), NodeArray(Nq,2), NodeArray(Nq,3), 'r.', 'MarkerSize', 20 )
%     daspect([1 1 1])
%     colorbar()
%     hold off
end