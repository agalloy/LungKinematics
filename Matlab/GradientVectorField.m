% Compute the gradient vector field on a surface from a scalar field.
% Inputs:
%   FaceArray = F X 3 interger array of face connectivity info
%   NodeArray = N X 3 float array of nodal/vertex positions
%   f = N X 1 float array giving the scalar functions value at each node
%   DEC = (Optional) A struct containing important geometric info about the
%       surface
% Outputs:
%   grad_f = N X 3 float array of gradient field vector components at each node 

function grad_f = GradientVectorField( FaceArray, NodeArray, f, DEC )
%% Assemble geometric prerequisites
    % Face and node numbers
    num_faces = size(FaceArray,1);

    % Determine if DEC needs to be computed or not
    if ~exist('DEC','var')
        DEC = AssembleDEC( FaceArray, NodeArray );
    end
    
    % Get info from DEC
    NodeAreas = DEC.NodeAreas;
    NodeFaceAreas = DEC.FaceNodeAreas';
    NodeFaceWeights = NodeFaceAreas ./ NodeAreas;
    
%% Compute the gradient field on each face and interpolate to nodes
    % Compute the potential gradient on every face
    grad_f = zeros(num_faces,3);
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
        ga = [ f(FaceArray(i,2)) - f(FaceArray(i,1))
               f(FaceArray(i,3)) - f(FaceArray(i,1)) ];
        % Get the gradient in the u1,u2 basis
        ga = CB*ga;
        % Get the gradient in the global basis
        grad_f(i,:) = (ga(1)*u1 + ga(2)*u2)';
    end
    
    % Interpolate face gradients onto each node
    grad_f = NodeFaceWeights * grad_f;
end