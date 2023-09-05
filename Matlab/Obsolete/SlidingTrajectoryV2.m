

function outStruct = SlidingTrajectoryV2(file,ref_lobe,invert_flag) 
%% Load in the FE mesh and displacement data
    % Get the subject from the file name
    [~, fname, ~] = fileparts(file);
    subject = extractBefore(fname,'_');
    
    % Load the mesh
    m_file = LoadMesh(file);
    load(m_file,'-mat','NodeArray','ElementArray','ZoneArray');
    
    % Load displacement data
    data = readmatrix(file);
    u = data(:,2:end);
    t_steps = size(u,2)/3;
    
    % Get the elements belonging to LLL and LUL
    el_start = ZoneArray(1,2) + 1;
    el_end = el_start + ZoneArray(2,2) - 1 ;
    el_LLL = ElementArray( el_start : el_end, : );
    el_start = el_end + 1;
    el_end = el_start + ZoneArray(3,2) - 1;
    el_LUL = ElementArray( el_start : el_end, : );   
    
    % Get only the surface nodes of each lobe
    fa_LLL = FESurface(el_LLL);
    ns_LLL = unique(fa_LLL);
    fa_LUL = FESurface(el_LUL);
    ns_LUL = unique(fa_LUL);
    
    if strcmp(ref_lobe,'LLL')
        % Set the reference and moving lobe
        ns_ref = ns_LLL;
        Xs_ref = NodeArray(ns_LLL,:);
        fa_ref = fa_LLL;
        ns_mov = ns_LUL;
        Xs_mov = NodeArray(ns_LUL,:);
        fa_mov = fa_LUL;
        
    elseif strcmp(ref_lobe,'LUL')
        % Set the reference and moving lobe
        ns_ref = ns_LUL;
        Xs_ref = NodeArray(ns_LUL,:);
        fa_ref = fa_LUL;
        ns_mov = ns_LLL;
        Xs_mov = NodeArray(ns_LLL,:);
        fa_mov = fa_LLL;
        
    else
        error('Invalid choice of reference lobe.')        
    end
    
%% Calculate sliding trajectory
% For each node on a give lobe's surface, find the nearest node
% on the other lobe's surface (within a given tolerance)
    xs_mov = Xs_mov;
    xs_ref = Xs_ref;
    
    for i = 2:t_steps
        t = (i-1)*3 + 1;
        
    % Locate pairs of nodes after one time step has passed and aligned the 
    % lobar surfaces better.
        if i == 3
        % Find closest points on the reference lobe surface to each node 
        % on the moving lobe surface (in the current configuration)
            cp_struct = ClosestPointTriSurf(fa_ref,NodeArray,xs_mov);
            
        % Get only the nodes near the reference surface or inside it
            tol = 0.5;
            dist = cp_struct.dist;
            inside = cp_struct.inside;
            pair = dist < tol | inside;
            num_pair = sum(pair);

            u_inc = zeros(num_pair,t_steps*3);
            
        % Get the nodal coordinates of paired points and the cooresponding
        % face and interpolation weights on the reference surface
            np_mov = ns_mov(pair);
            Xp_mov = Xs_mov(pair,:);
            face_cp = cp_struct.face(pair);
            Xcp_ref = cp_struct.Xcp(pair,:);
            Ncp = cp_struct.Ncp(pair,:);
            
        % Initialize the displacement array for reference surface points
            up_ref = zeros(num_pair,3);
        end
    
    % Calculate the incremental displacements for each node pair identified
    % in the first deformed time step
        if i >= 3
        % Find the incremental displacement for the moving lobe points
            up_mov = u(np_mov,t:t+2) - u(np_mov,t-3:t-1);
        
        % Find the incremental displacement for the reference lobe points
            for j = 1:num_pair
                % Get the nodes on the closest point's face and their
                % incremental displacements
                nf = fa_ref(face_cp(j),:);
                up_nf = u(nf,t:t+2) - u(nf,t-3:t-1);
                
                % Interpolate the incremental displacement at the closest
                % point
                up_ref(j,:) = Ncp(j,:)*up_nf;
            end

        % Get the relative incremental displacement betweeh the paired points
            u_inc(:,t:t+2) = up_mov - up_ref;
        end
        
    % Update the current configuration
        xs_mov = xs_mov + u(ns_mov,t:t+2) - u(ns_mov,t-3:t-1);
        xs_ref = xs_ref + u(ns_ref,t:t+2) - u(ns_ref,t-3:t-1);
    end
    
    % Get the incremental arc length for each trajectory
    s_inc = zeros(num_pair,t_steps);
    u_total = 0;
    traj = repmat(Xp_mov,1,t_steps);
    for i = 1:t_steps
        t = (i-1)*3 + 1;
        s_inc(:,i) = vecnorm(u_inc(:,t:t+2),2,2);
        u_total = u_total + u_inc(:,t:t+2);
        traj(:,t:t+2) = traj(:,t:t+2) + u_total;
    end
    
    % Add up the incremental arc lengths for the total arc length
    s_total = sum(s_inc,2);
    
    % Get the "straight line" arc-length
    s_line = sqrt(sum(u_total.^2,2));
    
%% Generate Output Structure

    outStruct.NodeArray = NodeArray;
    outStruct.ElementArray = ElementArray;
    outStruct.ZoneArray = ZoneArray;
    outStruct.nodes_mov = np_mov;
    outStruct.s_inc = s_inc;
    outStruct.s_total = s_total;
    outStruct.s_line = s_line;
    outStruct.traj = traj;
        
%% Plot Sliding Trajectories Against Moving Lobe Mesh
    % To reduce figure clutter sample to a given spacing
    n_samp = GridSample(Xp_mov,15);
    %n_samp = 1:length(Xp_mov);
    
    % Plot sliding trajectories against mesh
    figure()
    hold on
    for i = 1:length(n_samp)
        a = traj(n_samp(i),:);
        a = reshape(a,3,[]);
        a = a';
        plot3(a(:,1),a(:,2),a(:,3),'LineWidth',2)
    end
    plot3(Xp_mov(n_samp,1),Xp_mov(n_samp,2),Xp_mov(n_samp,3),'.k','MarkerSize',15)
    trimesh(fa_mov,NodeArray(:,1),NodeArray(:,2),NodeArray(:,3), ...
        'EdgeColor',[0.5,0.5,0.5])
    title(sprintf('Sliding Trajectories: %s',subject))
    
    % Camera control
    daspect([1 1 1])
    if invert_flag
        set(gca, 'Zdir', 'reverse')
        set(gca, 'Ydir', 'reverse')
    end
    
    hold off
end