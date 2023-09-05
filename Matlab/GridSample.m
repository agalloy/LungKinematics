% Return a set of query points sampled (roughly) uniformly with a given spacing

function n_resamp = GridSample(Xq,spacing)
    % Get the bounds of the query points
    bbox_min = min(Xq,[],1);
    bbox_max = max(Xq,[],1);
    
    % Create a set of grid points within these bounds with the given
    % spacing
    Xg_x = bbox_min(1): spacing: bbox_max(1);
    Xg_y = bbox_min(2): spacing: bbox_max(2);
    Xg_z = bbox_min(3): spacing: bbox_max(3);
    [Xg_x, Xg_y, Xg_z] = meshgrid(Xg_x,Xg_y,Xg_z);
    Xg = [ reshape(Xg_x,[],1), reshape(Xg_y,[],1), reshape(Xg_z,[],1)];
    
    % Get the closest query point to each grid point
    [k, dist] = dsearchn(Xq,Xg);
    
    % Find the grid points that have a distance of less than the
    % spacing to a query point and return the indices to those query points
    close = dist < spacing;
    n_resamp = k(close);    
end