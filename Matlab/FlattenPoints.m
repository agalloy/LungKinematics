% Flatten a point cloud along its axis of least variation

function V_flat = FlattenPoints( V )
    % Shift V so that the origin is at the mean
    V_shift = V - mean(V,1);
    % Perform singular value decomposition to get principal axes
    [~,~,s_axes] = svd( V_shift );
   
    % Rewrite V in terms of its principal axes
    V_flat = V_shift * s_axes;
    % Remove the component of V along the minimum axis
    V_flat(:,3) = 0;
end