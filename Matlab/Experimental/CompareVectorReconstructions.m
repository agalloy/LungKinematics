% Compare different methods of reconstructing a vector field
% Inputs:
%   FaceArray = F X 3 interger array of face connectivity info
%   NodeArray = N X 3 float array of nodal/vertex positions
%   X = N X 3 float array of nodal vector field components
%   DEC = (Optional) A struct containing important geometric info about the
%       surface
function hlink = CompareVectorReconstructions( FaceArray, NodeArray, X, DEC )
%% Assemble geometric prerequisites
    if all(NodeArray(:,3) == 0)
        ndim = 2;
    else
        ndim = 3;
    end
    % Determine if DEC needs to be computed or not
    if ~exist('DEC','var')
        DEC = AssembleDEC( FaceArray, NodeArray );
    end
    
    % Get info from DEC
    EdgeArray = DEC.EdgeArray;
    EdgeLengths = DEC.EdgeLengths;
    EdgeDir = DEC.EdgeDir;
    EdgeVec = EdgeDir .* EdgeLengths;
    NodeEdgeAreas = DEC.NodeEdgeAreas;
    hs0 = DEC.hs0;
    hs1 = DEC.hs1;
    
    % Face and node numbers
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
    
    X_tall = reshape(X,[],1);
    omega = VecToForm * X_tall;

%% Reconstruct the original vector field from the 1-form
    % Convert that 1-form to the vector componets along edges
    Y_edge = omega .* EdgeDir ./ EdgeLengths;
    
    % Interpolate the edge vector components to node vectors
    NodeEdgeWeights = NodeEdgeAreas ./ sum(NodeEdgeAreas,2);
    Y = NodeEdgeWeights * Y_edge;
        
%% Reconstruct the original vector field from the scalar field
    HHD_struct = HHD_GradientRecon( FaceArray, NodeArray, X );
    Z = HHD_struct.diff_alpha + HHD_struct.codiff_beta;
    
%% Assemble a linear operator that turns a 1-form to a vector field
    % Use the pseudo-inverse to be sure the vector field is as similar to
    % original as possible
    W_tall = lsqminnorm( VecToForm, omega, 'warn' );
    W = reshape( W_tall, [], 3 );
    
%% Compare vector fields reconstructed from both methods visually
    options = struct();
    
    figure()
    sp_array(1) = subplot(2,2,1);
    title('Original Vector Field')
    PlotTriSurfStreamline( FaceArray, NodeArray(:,1:ndim), X(:,1:ndim), options )
    daspect([1 1 1])
    %figure()
    sp_array(2) = subplot(2,2,2);
    title('Node interpolated 1-form reconstruction')
    PlotTriSurfStreamline( FaceArray, NodeArray(:,1:ndim), Y(:,1:ndim), options )
    daspect([1 1 1])
    sp_array(3) = subplot(2,2,3);
    title('Least squares 1-form reconstruction')
    PlotTriSurfStreamline( FaceArray, NodeArray(:,1:ndim), Y(:,1:ndim), options )
    daspect([1 1 1])
    %figure()
    sp_array(4) = subplot(2,2,4);
    title('Gradient reconstruction')
    PlotTriSurfStreamline( FaceArray, NodeArray(:,1:ndim), Z(:,1:ndim), options )
    daspect([1 1 1])
    
    hlink = linkprop( sp_array, {'CameraPosition','CameraUpVector'} );
    
%% Compare vector fields quantitatively
    X_reshape = reshape(X(:,1:ndim),[],1);
    Y_reshape = reshape(Y(:,1:ndim),[],1);
    Z_reshape = reshape(Z(:,1:ndim),[],1);
    W_reshape = reshape(W(:,1:ndim),[],1);

    
    residual_Y = abs( Y_reshape - X_reshape ) ./ abs( X_reshape );
    residual_Z = abs( Z_reshape - X_reshape ) ./ abs( X_reshape );
    residual_W = abs( W_reshape - X_reshape ) ./ abs( X_reshape );
    disp( 'Vector field reconstruction comparison:' )
    fprintf( '\n Max Node Interpolated Recon normalized residual: %f\n', max(residual_Y) )
    fprintf( '\n Max Least Squares Recon normalized residual: %f\n', max(residual_W) )
    fprintf( '\n Max Gradient Recon normalized residual: %f\n', max(residual_Z) )
    fprintf( '\n Mean Node Interpolated Recon normalized residual: %f\n', mean(residual_Y) )
    fprintf( '\n Mean Least Squares Recon normalized residual: %f\n', mean(residual_W) )
    fprintf( '\n Mean Gradient Recon normalized residual: %f\n', mean(residual_Z) )
    fprintf( '\n' )
    
%% Compare vector field and 1-form L2-Norms (magnitdues)
    % Inner product operator on nodal vector fields
    IPv = hs0;    
    % Inner product operator on 1-forms
    IP1 = hs1;
    
    % Compute L2 Norms
    L2_X = sqrt(sum(dot( X, IPv*X )));
    L2_Y = sqrt(sum(dot( Y, IPv*Y )));
    L2_Z = sqrt(sum(dot( Z, IPv*Z )));
    L2_W = sqrt(sum(dot( W, IPv*W )));
    L2_omega = sqrt( omega' * IP1 * omega );
    
    L2_table = table( L2_X, L2_omega, L2_Y, L2_W, L2_Z,...
                      'VariableNames',["Original","Intermediate 1-form", "Node Interpolated Recon", "Least Squares Recon", "Gradient Recon"] );
    disp('Vector field L2-norms')
    disp(L2_table)
end