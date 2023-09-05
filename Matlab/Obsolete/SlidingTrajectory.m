

function outStruct = SlidingTrajectory(file,ref_lobe,invert_flag) 
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
        ns_mov = ns_LUL;
        Xs_mov = NodeArray(ns_LUL,:);
        fa_mov = fa_LUL;
        
    elseif strcmp(ref_lobe,'LUL')
        % Set the reference and moving lobe
        ns_ref = ns_LUL;
        Xs_ref = NodeArray(ns_LUL,:);
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
        % Find closest node on the reference lobe to each node on the moving lobe
        % (in the current configuration)
            k = dsearchn(xs_ref,xs_mov);

        % Remove node pairs whose distance does not meet a tolerance
            tol = 1;
            p_dist = xs_mov - xs_ref(k,:);
            p_dist = vecnorm(p_dist,2,2);
            pair = p_dist < tol;
            num_pair = sum(pair);

            u_inc = zeros(num_pair,t_steps*3);
        
        %Get nodal coordinates of paired nodes
            np_mov = ns_mov(pair);
            Xp_mov = Xs_mov(pair,:);
            np_ref = ns_ref(k(pair));
            Xp_ref = Xs_ref(k(pair),:);
        end
    
    % Calculate the incremental displacements for each node pair identified
    % in the first deformed time step
        if i >= 3
        % Find the incremental displacement vectors for each node in a pair 
        % then find the difference in the incremental displacement
            up_mov = u(np_mov,t:t+2) - u(np_mov,t-3:t-1);
            up_ref = u(np_ref,t:t+2) - u(np_ref,t-3:t-1);    
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
    outStruct.nodes_ref = np_ref;
    outStruct.nodes_mov = np_mov;
    outStruct.s_inc = s_inc;
    outStruct.s_total = s_total;
    outStruct.s_line = s_line;
    outStruct.traj = traj;
        
%% Plot Sliding Trajectories Against Moving Lobe Mesh
    
    % Plot sliding trajectories against mesh
    figure()
    hold on
    for i = 1:num_pair
        a = traj(i,:);
        a = reshape(a,3,[]);
        a = a';
        plot3(a(:,1),a(:,2),a(:,3),'LineWidth',2)
    end
    plot3(Xp_mov(:,1),Xp_mov(:,2),Xp_mov(:,3),'.k','MarkerSize',15)
    trimesh(fa_mov,NodeArray(:,1),NodeArray(:,2),NodeArray(:,3), ...
        'EdgeColor',[0.5,0.5,0.5],'FaceAlpha',0.5)
    title(sprintf('Sliding Trajectories: %s',subject))
    
    % Camera control
    daspect([1 1 1])
    if invert_flag
        set(gca, 'Zdir', 'reverse')
        set(gca, 'Ydir', 'reverse')
    end
    
    hold off
end