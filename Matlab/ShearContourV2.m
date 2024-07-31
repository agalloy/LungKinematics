% Interpolate FE mesh and displacements onto a voxel grid with elements
% that cross between solid domains and then calculate max shear on new grid.


function [avg_maxshear, out_struct] = ShearContourV2( NodeArray, ElementArray, eID, DispArray, spacing )
%% Step 1: Interpolate FE Results onto a voxel grid
    % Create array of points to sample
    s_start = (floor(min(NodeArray)./spacing)-1).*spacing;
    s_end = (ceil(max(NodeArray)./spacing)+1).*spacing;
    x = s_start(1):spacing(1):s_end(1);
    y = s_start(2):spacing(2):s_end(2);
    z = s_start(3):spacing(3):s_end(3);
    grid_size = [ size(x,2), size(y,2), size(z,2) ];
    [grid_x, grid_y, grid_z] = meshgrid(x,y,z); 
    grid_points = [ reshape(grid_x,[],1), reshape(grid_y,[],1), reshape(grid_z,[],1) ];
    [grid_i, grid_j, grid_k] = meshgrid( 1:grid_size(1), 1:grid_size(2), 1:grid_size(3) );
    grid_index = sub2ind( grid_size, reshape(grid_i,[],1), reshape(grid_j,[],1), reshape(grid_k,[],1) );
    
    % Interpolate displacements to a voxel grid
    DomainMask = zeros(grid_size); % Image mask of mesh domains
    DispXImage = nan(grid_size); % X-Displacements interpolated to new grid
    DispYImage = nan(grid_size); % Y-Displacements interpolated to new grid
    DispZImage = nan(grid_size); % Z-Displacements interpolated to new grid

    % Loop through each solid domain
    warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
    for i = 2:max(eID)
        % Convert mesh into a triangulation object
        triang = triangulation( ElementArray(eID==i,:), NodeArray );
        
        % Figure out if each sample point is inside the triangulation
        % If so provide an an element ID and shape functions 
        [ID,N] = pointLocation(triang,grid_points);
        in_elem = ~isnan(ID);
        
        % Set the mask to appropriate value if a sample point is inside
        % the domain
        DomainMask( grid_index(in_elem) ) = i;

        % Use the shape functions to interpolate displacement values
        domain_EA = ElementArray( eID==i, : );
        nodal_dispX = DispArray( domain_EA(ID(in_elem),:), 1 );
        nodal_dispX = reshape(nodal_dispX,[],4);
        DispXImage( grid_index(in_elem) ) = dot( N(in_elem,:), nodal_dispX, 2 );

        nodal_dispY = DispArray( domain_EA(ID(in_elem),:), 2 );
        nodal_dispY = reshape(nodal_dispY,[],4);
        DispYImage( grid_index(in_elem) ) = dot( N(in_elem,:), nodal_dispY, 2 );

        nodal_dispZ = DispArray( domain_EA(ID(in_elem),:), 3 );
        nodal_dispZ = reshape(nodal_dispZ,[],4);
        DispZImage( grid_index(in_elem) ) = dot( N(in_elem,:), nodal_dispZ, 2 );
    end
    warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
    
%% Step 2: Crete a new brick8 mesh from the new voxel image
    % March through each grid point and check if it can be part of an element.
    % If so, add the node to the new node array (NodeArray2) checking for
    % redundancies as well as adding the element to the new element array
    
    % Pre-allocate new mesh arrays
    num_points = prod(grid_size);
    NodeCatalog = zeros(grid_size);
    NodeArray2 = zeros(num_points,3);
    DispArray2 = zeros(num_points,3);
    ElementArray2 = zeros(num_points,8);
    FissureIndex = zeros(num_points,1);
    
    % Get only object voxels to speed things up
    indices = find( DomainMask );
    [I,J,K] = ind2sub( grid_size, indices );
    
    % Loop through every voxel in mask and check if an element can be made
    % from it
    count_n = 0; %initialize node count
    count_e = 0; %initialize element count
    % Set up neighborhood relationships
    i_nbhd = zeros(2,2,2);
    i_nbhd(2,:,:) = 1;
    j_nbhd = zeros(2,2,2);
    j_nbhd(:,2,:) = 1;
    k_nbhd = zeros(2,2,2);
    k_nbhd(:,:,2) = 1;
    for n = 1:length(indices)
        i = I(n);
        j = J(n);
        k = K(n);
        %If all voxels in the brick are nonzero create new brick element
        if all( DomainMask(i:i+1,j:j+1,k:k+1) )
            % Find uncatalogged nodes in brick to add to Node Catalog
            NodeBrick = NodeCatalog(i:i+1,j:j+1,k:k+1);
            to_catalog = NodeBrick == 0;
            num_new_nodes = nnz(to_catalog);
            new_nodes = (count_n + 1) : (count_n + num_new_nodes);
            count_n = count_n + num_new_nodes;
            NodeBrick(to_catalog) = new_nodes;
            NodeCatalog(i:i+1,j:j+1,k:k+1) = NodeBrick;
            
            % Update NodeArray2
            i_new = i + i_nbhd(to_catalog);
            j_new = j + j_nbhd(to_catalog);
            k_new = k + k_nbhd(to_catalog);
            NodeArray2(new_nodes,:) = [ x(i_new)', y(j_new)', z(k_new)' ];
            
            % Update DispArray2
            DispXBrick = DispXImage(i:i+1,j:j+1,k:k+1);
            DispYBrick = DispYImage(i:i+1,j:j+1,k:k+1);
            DispZBrick = DispZImage(i:i+1,j:j+1,k:k+1);
            DispArray2(new_nodes,:) = [ DispXBrick(to_catalog), DispYBrick(to_catalog), DispZBrick(to_catalog) ];
                
            % Update ElementArray2
            count_e = count_e + 1;
            ElementArray2(count_e,:) = reshape( NodeBrick, [], 1 );
            
            % Update FissureIndex
            unique_domains = unique( DomainMask(i:i+1,j:j+1,k:k+1) );
            if numel(unique_domains) == 1
                FissureIndex(count_e) = 0;
            else
                FissureIndex(count_e) = 1;
            end
        end
    end
    
    %Clip the unneeded portions of NodeArray2, ElementArray2, and
    %FissureIndex
    NodeArray2 = NodeArray2(1:count_n,:);
    DispArray2 = DispArray2(1:count_n,:);
    ElementArray2 = ElementArray2(1:count_e,:);
    FissureIndex = FissureIndex(1:count_e);
%% Step 3: Calculate Max Shear at each element of the new mesh
    %Create a matrix with the derivatives of the brick8 shape functions
    nodes = ElementArray2(1,:);
%     dNdS = 1/4*[
%         -1 -1 -1
%          1 -1 -1
%          1  1 -1
%         -1  1 -1
%         -1 -1  1
%          1 -1  1
%          1  1  1
%         -1  1  1
%     ];
    dNdS = 1/4*[
        -1 -1 -1
         1 -1 -1
        -1  1 -1
         1  1 -1
        -1 -1  1
         1 -1  1
        -1  1  1
         1  1  1
    ];
    dXdS = NodeArray2(nodes,:)'*dNdS;
    dNidXj = dNdS*(dXdS)^-1;
    
    %Calculate the max shear for each element
    MaxShear = zeros(size(ElementArray2,1),1);
    for i = 1:size(ElementArray2,1)
        %Generate a matrix of each node of the element's displacements
        u = [DispArray2(ElementArray2(i,1),:)',DispArray2(ElementArray2(i,2),:)',DispArray2(ElementArray2(i,3),:)',DispArray2(ElementArray2(i,4),:)',DispArray2(ElementArray2(i,5),:)',DispArray2(ElementArray2(i,6),:)',DispArray2(ElementArray2(i,7),:)',DispArray2(ElementArray2(i,8),:)'];
        %Calculate deformation gradient dxi/dXj
        F = u*dNidXj+eye(3);
        %Calculate Right Cauchy Green Tensor
        C = F'*F;
        %Principal Stretches from C
        lambda_sq = sort(eig(C), 'descend');
        lambda = sqrt(lambda_sq);
        %Max shear from principal stretches
        MaxShear(i) = (lambda(1)-lambda(3))/2;
    end
%% Assemble outputs
    % Calculate average max shear on fissures
    avg_maxshear = mean( MaxShear(FissureIndex == 1) );
    
    % Additional output structure with all info
    out_struct = struct();
    out_struct.NodeArray = NodeArray2;
    out_struct.ElementArray = ElementArray2;
    out_struct.FissureIndex = FissureIndex;
    out_struct.MaxShear = MaxShear;
end