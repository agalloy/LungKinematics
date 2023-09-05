% Get the closest points on a triangulated surface to the query points.

function outStruct = ClosestPointTriSurf( FaceArray, NodeArray, Xq )
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
    
    % Get the surface nodes and their positions
    ns = unique( FaceArray );
    Xs = NodeArray(ns,:);
    
    % Get the node on the surface closest to the query point
    k = dsearchn(Xs,Xq);
    
%% Loop through each face to find its distance to each query point
    % Set temporary values for variables before finding closest point
    % values.
    Xcp = zeros(length(Xq),3); % Coordinates of cp
    dist = zeros(length(Xq),1) + 1000000; % Distance between query point and cp
    face = zeros(length(Xq),1);

    Xcp_temp = zeros(length(Xq),3) + 10000000;
    dist_temp = zeros(length(Xq),1);
    
    % Loop through each face on the surface
    for i = 1:length(FaceArray)        
        % Project query nodes onto the same plane as the current face
        Xq_rel = Xq - NodeArray(FaceArray(i,1),:);
        Xq_pro = Xq_rel - (Xq_rel*normals(i,:)')*normals(i,:) + NodeArray(FaceArray(i,1),:);
        
        % Check if the projected point is in the triangle
        in12 = ( Xq_pro - NodeArray(FaceArray(i,1),:) ) * norm12(i,:)' > 0;
        in23 = ( Xq_pro - NodeArray(FaceArray(i,2),:) ) * norm23(i,:)' > 0;
        in31 = ( Xq_pro - NodeArray(FaceArray(i,3),:) ) * norm31(i,:)' > 0;
        intri = all( [in12, in23, in31], 2 );
        
        % Get closest point between each edge and the projected point
        cp12 = ( Xq_pro - NodeArray(FaceArray(i,1),:) ) * edge12(i,:)' ./ norm(edge12(i,:),2)^2;
        cp12(cp12 > 1) = 1;
        cp12(cp12 < 0) = 0;
        cp12 = cp12*edge12(i,:) + NodeArray(FaceArray(i,1),:);
        
        cp23 = ( Xq_pro - NodeArray(FaceArray(i,2),:) ) * edge23(i,:)' ./ norm(edge23(i,:),2)^2;
        cp23(cp23 > 1) = 1;
        cp23(cp23 < 0) = 0;
        cp23 = cp23*edge23(i,:) + NodeArray(FaceArray(i,2),:);
        
        cp31 = ( Xq_pro - NodeArray(FaceArray(i,3),:) ) * edge31(i,:)' ./ norm(edge31(i,:),2)^2;
        cp31(cp31 > 1) = 1;
        cp31(cp31 < 0) = 0;
        cp31 = cp31*edge31(i,:) + NodeArray(FaceArray(i,3),:);
        
        % Find the distance between these edge points and the query point
        dcp12 = vecnorm( Xq - cp12, 2, 2);
        dcp23 = vecnorm( Xq - cp23, 2, 2);
        dcp31 = vecnorm( Xq - cp31, 2, 2);
        
        % If the projected point is inside the face, then that is the closest
        % point on the face to Xq.
        Xcp_temp(intri,:) = Xq_pro(intri,:);
        dist_temp(intri) = vecnorm( Xq(intri,:) - Xcp_temp(intri,:), 2, 2);
        
        % If the projected point is outside of the face, then the closest
        % point is one of the edge points above
        edge_array = [ dcp12(~intri), dcp23(~intri), dcp31(~intri)];
        edgep_array = [ cp12(~intri,:), cp23(~intri,:), cp31(~intri,:)];
        [dist_temp(~intri), edgeind] = min(edge_array,[],2);
        Xcp_temp(~intri,:) = edgep_array( :, edgeind*3-2 : edgeind*3);
        
        % If this face is the closest point detected for a query node so far
        % save this distance.
        Xcp(dist_temp < dist,:) = Xcp_temp(dist_temp < dist,:);
        face(dist_temp < dist) = i;
        dist(dist_temp < dist) = dist_temp(dist_temp < dist);
                
    end
    
%% Loop through each node to calculate shape functions and side
    Ncp = zeros(length(Xq),3);
    inside = false(length(Xq),1);
    
    for i = 1:length(Xq)
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