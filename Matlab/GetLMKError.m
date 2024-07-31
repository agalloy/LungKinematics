

function out_struct = GetLMKError(NodeArray, ElementArray, DispArray, lmk_file)
%% Read landmark data
    lmk_data = readmatrix(lmk_file);
    X_lmk = lmk_data(:,1:3);
    u_lmk = lmk_data(:,4:6);
    
%% Step 2: Interpolate shape functions to find the displacements at landmark positions
    %Create an array of nodal positions for each element
    TetArray = zeros(length(ElementArray),4,3);
    TetArray(:,1,:) = NodeArray(ElementArray(:,1),:);
    TetArray(:,2,:) = NodeArray(ElementArray(:,2),:);
    TetArray(:,3,:) = NodeArray(ElementArray(:,3),:);
    TetArray(:,4,:) = NodeArray(ElementArray(:,4),:);
    
    %Create an array of bounding box cooridinates for each element
    BoxArray = zeros(length(ElementArray),6);
    BoxArray(:,1) = min(TetArray(:,:,1),[],2);
    BoxArray(:,2) = max(TetArray(:,:,1),[],2);
    BoxArray(:,3) = min(TetArray(:,:,2),[],2);
    BoxArray(:,4) = max(TetArray(:,:,2),[],2);
    BoxArray(:,5) = min(TetArray(:,:,3),[],2);
    BoxArray(:,6) = max(TetArray(:,:,3),[],2);
    
    %Loop through each landmark and calculate Febio's predicted
    %displacement at the landmark position
    u_fea = zeros(size(X_lmk,1),3);
    for i = 1:size(X_lmk,1)
        %Get a list of elements the landmark is inside the bounding box of
        %to reduce the number of elements to run through
        in_box_x = BoxArray(:,1) <= X_lmk(i,1) & BoxArray(:,2) >= X_lmk(i,1);
        in_box_y = BoxArray(:,3) <= X_lmk(i,2) & BoxArray(:,4) >= X_lmk(i,2);
        in_box_z = BoxArray(:,5) <= X_lmk(i,3) & BoxArray(:,6) >= X_lmk(i,3);
        in_box = in_box_x & in_box_y & in_box_z;
        
        %Loop through each element in TetList
        for j = find(in_box)'
            %Create a linear map between the cartesian coordinates and
            %the element shape functions
            A = [
                1               1               1               1
                TetArray(j,1,1) TetArray(j,2,1) TetArray(j,3,1) TetArray(j,4,1)
                TetArray(j,1,2) TetArray(j,2,2) TetArray(j,3,2) TetArray(j,4,2)
                TetArray(j,1,3) TetArray(j,2,3) TetArray(j,3,3) TetArray(j,4,3)
            ];
            M = A^(-1);
            
            %Evaluate element shape functions at the landmark position
            N_lmk = M*[1,X_lmk(i,:)]';
            
            %If the landmark is inside the element, then 0 <= Ni <= 1 for
            %every node
            inside = N_lmk >= 0 & N_lmk <= 1;
            inside = all(inside);
            
            %If inside interpolate displacements to that point using shape
            %functions and find landmark error
            if inside
                if u_fea ~= 0
                    error('Duplicate!')
                end
                % Get element nodal displacements
                u_node = DispArray(ElementArray(j,:),:)';
                % Interpolate nodal displacements at landmark position
                u_fea(i,:) = (u_node*N_lmk)';
            end
        end
    end   


%% Assemble output structure
    out_struct.X_lmk = X_lmk;
    out_struct.u_lmk = u_lmk;
    out_struct.u_fea = u_fea;
end