% Get the closest points on a triangulated surface to the query points.

function outStruct = ClosestPointTriSurfV2( FaceArray, NodeArray, Xq )
%% Calculate some initial helpful variables
    % Get each face normal
    edge12 = NodeArray(FaceArray(:,2),:) - NodeArray(FaceArray(:,1),:);
    edge23 = NodeArray(FaceArray(:,3),:) - NodeArray(FaceArray(:,2),:);
    edge31 = NodeArray(FaceArray(:,1),:) - NodeArray(FaceArray(:,3),:);
    normals = cross( edge12, -edge31 );
    normals = normals ./ vecnorm(normals,2,2);
    
    % Get each edge normal
    norm12 = cross( normals, edge12 ) ./ vecnorm(cross( normals, edge12 ),2,2);
    norm23 = cross( normals, edge23 ) ./ vecnorm(cross( normals, edge23 ),2,2);
    norm31 = cross( normals, edge31 ) ./ vecnorm(cross( normals, edge31 ),2,2);
    
    % Get the maximum edge length of the surface
    l_edge12 = vecnorm(edge12,2,2);
    l_edge23 = vecnorm(edge23,2,2);
    l_edge31 = vecnorm(edge31,2,2);
    max_edge = max( [l_edge12; l_edge23; l_edge31] );
        
    % Get the surface nodes and their positions
    ns = unique( FaceArray );
    Xs = NodeArray(ns,:);
    
    % Get the node on the surface closest to the query point
    k = dsearchn(Xs,Xq);
    k = ns(k);
    X_k = NodeArray(k,:);
           
%% Loop through each query node to find the closest face
    % Set temporary values for variables before finding closest point
    % values.
    Xcp = nan(size(Xq,1),3); % Coordinates of cp
    dist = nan(size(Xq,1),1); % Distance between query point and cp
    face = zeros(size(Xq,1),1); % Nearest face to Xq
    Ncp = zeros(size(Xq,1),3); % Shape function values at Xcp
    inside = false(size(Xq,1),1); % Whether Xq is inside the surface
        
    % Loop through each query node
%    tic
    for i = 1:size(Xq,1)    
%         if ~mod(i,100)
%             fprintf('\n%i nodes remaining.\n',length(Xq) - i)
%             toc
%             tic
%         end
        
        % Rule out faces far from nearest node
        dist_k1 = vecnorm( NodeArray(FaceArray(:,1),:) - X_k(i,:), 2, 2 );
        near_face = find( dist_k1 <= 2*max_edge );
        
        % Get face vertices
        V1 = NodeArray( FaceArray(near_face,1), : );
        V2 = NodeArray( FaceArray(near_face,2), : );
        V3 = NodeArray( FaceArray(near_face,3), : );
        
        % Project query node onto the same plane as the faces
        Xq_rel = Xq(i,:) - V1;
        Xq_pro = Xq_rel - dot(Xq_rel,normals(near_face,:),2).*normals(near_face,:) + V1;
        
        % Check if the projected point is in the triangle
        in12 = dot( ( Xq_pro - V1 ), norm12(near_face,:), 2) > 0;
        in23 = dot( ( Xq_pro - V2 ), norm23(near_face,:), 2) > 0;
        in31 = dot( ( Xq_pro - V3 ), norm31(near_face,:), 2) > 0;
        intri = all( [in12, in23, in31], 2 );
        
        % Get closest point between each edge and the projected point
        cp12 = dot( ( Xq_pro - V1 ), edge12(near_face,:), 2) ./ vecnorm(edge12(near_face,:),2,2).^2;
        cp12(cp12 > 1) = 1;
        cp12(cp12 < 0) = 0;
        cp12 = cp12.*edge12(near_face,:) + V1;
        
        cp23 = dot( ( Xq_pro - V2 ), edge23(near_face,:), 2) ./ vecnorm(edge23(near_face,:),2,2).^2;
        cp23(cp23 > 1) = 1;
        cp23(cp23 < 0) = 0;
        cp23 = cp23.*edge23(near_face,:) + V2;
        
        cp31 = dot( ( Xq_pro - V3 ), edge31(near_face,:), 2) ./ vecnorm(edge31(near_face,:),2,2).^2;
        cp31(cp31 > 1) = 1;
        cp31(cp31 < 0) = 0;
        cp31 = cp31.*edge31(near_face,:) + V3;
        
        % Find the distance between these edge points and the query point
        dcp12 = vecnorm( Xq(i,:) - cp12, 2, 2);
        dcp23 = vecnorm( Xq(i,:) - cp23, 2, 2);
        dcp31 = vecnorm( Xq(i,:) - cp31, 2, 2);
        
        % If the projected point is inside the face, then that is the closest
        % point on the face to Xq.
        Xcp_temp = Xq_pro;
        dist_temp = vecnorm( Xq(i,:) - Xcp_temp, 2, 2);
        
        % If the projected point is outside of the face, then the closest
        % point is one of the edge points above
        edge_array = [ dcp12(~intri), dcp23(~intri), dcp31(~intri)];
        edgep_array = [ cp12(~intri,:), cp23(~intri,:), cp31(~intri,:)];
        [dist_temp(~intri), edgeind] = min(edge_array,[],2);
        Xcp_temp(~intri,:) = edgep_array( :, edgeind*3-2 : edgeind*3);
        
        % If this face is the closest point detected for a query node so far
        % save this distance.
        [dist(i), fID] = min( dist_temp );
        Xcp(i,:) = Xcp_temp( fID, : );
        face(i) = near_face(fID);
        
        % Calculate interpolation functions at the closest point projection
        f = face(i);
        M = [
            NodeArray(FaceArray(f,1),1) NodeArray(FaceArray(f,2),1) NodeArray(FaceArray(f,3),1)
            NodeArray(FaceArray(f,1),2) NodeArray(FaceArray(f,2),2) NodeArray(FaceArray(f,3),2)
            NodeArray(FaceArray(f,1),3) NodeArray(FaceArray(f,2),3) NodeArray(FaceArray(f,3),3)
            ]^-1;
        Ncp(i,:) = (M*Xcp(i,:)')';
        
        % Determine whether the query point is inside or outside the
        % surface using face normals
        inside(i,:) = (Xq(i,:) - NodeArray(FaceArray(f,1),:))*normals(f,:)' < 0;
    end

%% Generate errors
    if abs( sum(Ncp,2) - 1 ) > 10^-8 | any( Ncp < -10^-8, 'all')
        warning('Suspicious interpolation function values.')
    end

%% Assemble output
    % The distance between each query point and the reference surface
    outStruct.dist = dist;
    % The closest point projection of each query point onto the reference
    % surface
    outStruct.Xcp = Xcp;
    % The face on the reference surface closest to each query point
    outStruct.face = face;
    % The interpolation function weights at each closest point projection 
    outStruct.Ncp = Ncp;
    % The normal vecctors for each face on the reference surface
    outStruct.normals = normals;
    % Boolean indicator of wether a query point is inside or outside the
    % reference surface
    outStruct.inside = inside;
end