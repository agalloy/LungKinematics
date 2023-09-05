

function outStruct = RayTracing( FaceArray, NodeArray, Xq, norm_q )
%% Calculate some initial helpful variables
    % Get each face's vertex coordinates
    V1 = NodeArray(FaceArray(:,1),:);
    V2 = NodeArray(FaceArray(:,2),:);
    V3 = NodeArray(FaceArray(:,3),:);
    % Get each face edge vector
    edge12 = V2 - V1;
    edge23 = V3 - V2;
    edge31 = V1 - V3;
    edge12n = edge12 ./ vecnorm(edge12,2,2);
    edge23n = edge23 ./ vecnorm(edge23,2,2);
    edge31n = edge31 ./ vecnorm(edge31,2,2);
    % Get each reference face normal vector
    norm_ref = cross( edge12, -edge31 );
    norm_ref = norm_ref ./ vecnorm(norm_ref,2,2);    
    % Get each edge normal vector (in the face's plane)
    norm12 = cross( norm_ref, edge12 ) ./ vecnorm(cross( norm_ref, edge12 ),2,2);
    norm23 = cross( norm_ref, edge23 ) ./ vecnorm(cross( norm_ref, edge23 ),2,2);
    norm31 = cross( norm_ref, edge31 ) ./ vecnorm(cross( norm_ref, edge31 ),2,2);
    
%% Loop through each query point and find projection point
    % Initialize outputs
    face = nan(size(Xq,1),1);
    dist = nan(size(Xq,1),1);
    Xcp = nan(size(Xq,1),3);
    Ncp = nan(size(Xq,1),3);
    inside = false(size(Xq,1),1);
    
    % Loop trough each query point
    for i = 1:size(Xq,1)
    % Find the intersection between the normal line extending from the 
    % query point and each reference facet
        % Get the normal distance of Xq from the facet planes perspective
        Xq_rel = Xq(i,:) - V1;
        Xq_ndist_ref = dot(norm_ref,Xq_rel,2);
        % Get the distance of the facet plane along Xq's normal line
        dist_temp = Xq_ndist_ref ./ (norm_ref*norm_q(i,:)');
        % Get the intersection point along the normal line 
        Xcp_temp = Xq(i,:) - dist_temp.*norm_q(i,:);
        
    % Determine if the intersection point is within the triangle
        in12 = dot( ( Xcp_temp - V1 ), norm12, 2) > 0;
        in23 = dot( ( Xcp_temp - V2 ), norm23, 2) > 0;
        in31 = dot( ( Xcp_temp - V3 ), norm31, 2) > 0;
        intri = find( all( [in12, in23, in31], 2 ) );
        
    % Store projected points and distances as necessary
        if any(intri)
            % Store the values of the closest projection point
            [dist(i),a] = min(abs(dist_temp(intri)));
            Xcp(i,:) = Xcp_temp(intri(a),:);
            face(i) = intri(a);
            inside(i) = dist_temp(intri(a)) > 0;
            
            % Get the shape functions at the closest projection point
            f = face(i);
            M = [
                NodeArray(FaceArray(f,1),1) NodeArray(FaceArray(f,2),1) NodeArray(FaceArray(f,3),1)
                NodeArray(FaceArray(f,1),2) NodeArray(FaceArray(f,2),2) NodeArray(FaceArray(f,3),2)
                NodeArray(FaceArray(f,1),3) NodeArray(FaceArray(f,2),3) NodeArray(FaceArray(f,3),3)
                ]^-1;
            Ncp(i,:) = (M*Xcp(i,:)')';
        end
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
    outStruct.normals = norm_q;
    % Boolean indicator of wether a query point is inside or outside the
    % reference surface
    outStruct.inside = inside;
end