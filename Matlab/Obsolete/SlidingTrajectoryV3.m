

function outStruct = SlidingTrajectoryV3(inStruct) 
%% Load in inStruct parameters
    if isfield(inStruct,'dist_tol')
        dist_tol = inStruct.dist_tol;
    else
        dist_tol = 0.5;
    end    
    if isfield(inStruct,'plot_hide')
        plot_hide = inStruct.plot_hide;
    else
        plot_hide = false;
    end    
    if isfield(inStruct,'plot_scale')
        plot_scale = inStruct.plot_scale;
    else
        plot_scale = 1;
    end    
    if isfield(inStruct,'plot_invert')
        plot_invert = inStruct.plot_invert;
    else
        plot_invert = false;
    end    
    if isfield(inStruct,'plot_blank')
        plot_blank = inStruct.plot_blank;
    else
        plot_blank = -1;
    end    
    if isfield(inStruct,'end_step')
        end_step = inStruct.end_step;
    else
        end_step = 0;
    end    
    if isfield(inStruct,'start_step')
        start_step = inStruct.start_step;
    else
        start_step = 1;
    end    
    if isfield(inStruct,'two_step')
        two_step = inStruct.two_step;
    else
        two_step = false;
    end

%% Load in the FE mesh and displacement data
    if isfield(inStruct,'feb_file')
        % Load .feb file string
        feb_file = inStruct.feb_file;
        fID = fopen(feb_file);
        feb_txt = fscanf(fID,'%c');
        fclose(fID);

        % Get element sections (solid elements only) from .feb text
        elem_txt = extractBetween(feb_txt,'<Elements type="tet4"', ...
                                 sprintf('\r\n\t\t</Elements>') );

        % Loop through each element domain and fill up Element connectivity array
        ElementArray = cell( size(elem_txt,1), 1 );
        for i = 1:length(elem_txt)
            % Split text into rows and remove junk
            elem_txt_split = split(elem_txt{i},sprintf('\r\n\t\t\t'));
            elem_txt_split = elem_txt_split(2:end);
            elem_txt_split = strjoin(elem_txt_split);

            % Parse element connectivity from text
            EA_temp = sscanf(char(elem_txt_split),'<elem id="%d">%d,%d,%d,%d</elem> ');
            EA_temp = ( reshape(EA_temp,5,[]) )';
            EA_temp = EA_temp(:,2:end);    
            ElementArray{i} = EA_temp;
        end
        clear EA_temp elem_txt
        
        % Get node sections from the .feb text
        node_txt = extractBetween(feb_txt,'<Nodes', ...
                                 sprintf('\r\n\t\t</Nodes>') );
                             
        % Loop through each node domain and fill up Node position array
        NodeArray = cell( size(node_txt,1), 1 );
        for i = 1:length(node_txt)
            % Split text into rows and remove junk
            node_txt_split = split(node_txt{i},sprintf('\r\n\t\t\t'));
            node_txt_split = node_txt_split(2:end);
            node_txt_split = strjoin(node_txt_split);

            % Parse node positions from text
            NA_temp = sscanf(char(node_txt_split),'<node id="%d">%f,%f,%f</node> ');
            NA_temp = ( reshape(NA_temp,4,[]) )';
            NA_temp = NA_temp(:,2:end);    
            NodeArray{i} = NA_temp;
        end
        NodeArray = cell2mat(NodeArray);
        clear NA_temp node_txt feb_txt
        
        if isfield(inStruct,'pos_data')
            % Read FEA position data
            pos_file = inStruct.pos_data;
            pos_data = readmatrix(pos_file);
            pos_data = pos_data(:,2:end);
            t_steps = size(pos_data,2) / 3;

            % Set the appropriate reference config based on the start step
            NodeArray = pos_data( :, (start_step*3-2):(start_step*3) );
            % Get the displacement cell array for all requested time steps
            u = pos_data - repmat(NodeArray,1,t_steps);
            u = mat2cell( u, size(u,1), repmat(3,1,t_steps) );
            if end_step ~= 0
                u = u(start_step:end_step);
                t_steps = end_step - start_step + 1;
            end                  
        elseif isfield(inStruct,'disp_data')
            % Read FEA disp data
            disp_file = inStruct.disp_data;
            disp_data = readmatrix(disp_file);
            disp_data = disp_data(:,2:end);
            t_steps = size(disp_data,2) / 3;

            % Set the appropriate reference config based on start_step
            u_offset = disp_data( :, (start_step*3-2):(start_step*3) );
            NodeArray = NodeArray + u_offset;
            u = disp_data - repmat(u_offset,1,t_steps);
            u = mat2cell( u, size(u,1), repmat(3,1,t_steps) );
            % Get the displacement cell array for all requested time steps
            if end_step ~= 0
                u = u(start_step:end_step);
                t_steps = end_step - start_step + 1;
            end      
        else
            error('Position or displacement data required!')
        end
        
        % If 2-step is set, only consider the first and final config
        if two_step
            u = [ u(1); u(end)];
            t_steps = 2;
        end
        
        % Set the reference and moving surfaces
        ns_ref = unique( cell2mat(ElementArray(inStruct.ref_surface)) );
        Xs_ref = NodeArray(ns_ref,:);
        fa_ref = FESurface( cell2mat(ElementArray(inStruct.ref_surface)) );
        ns_mov = unique( cell2mat(ElementArray(inStruct.mov_surface)) );
        Xs_mov = NodeArray(ns_mov,:);
        fa_mov = FESurface( cell2mat(ElementArray(inStruct.mov_surface)) );
        ElementArray = cell2mat(ElementArray);        
    else
        error('FEA Mesh must be provided to obtain a sliding trajectory.')
    end
    
%% Calculate sliding trajectory
    num_mov_nodes = size(ns_mov,1);
    num_ref_nodes = size(ns_ref,1);
    
    % Initialize the moving point spatial coordinate array
    ut_inc = nan(size(NodeArray,1),3*t_steps - 3);
    ut_inc = mat2cell( ut_inc, size(NodeArray,1), repmat(3,1,t_steps-1) );
    % Initialize the closest reference surface point (at the given time 
    % step) reference coordinate array
    Xcp = nan(size(NodeArray,1),3*t_steps);
    Xcp = mat2cell( Xcp, size(NodeArray,1), repmat(3,1,t_steps) );
    
    % Initialize the closest reference surface point displacement array
    Ucp = nan(size(NodeArray,1),3*t_steps);
    Ucp = mat2cell( Ucp, size(NodeArray,1), repmat(3,1,t_steps) );
    
    % Initialize is_sliding list
    is_sliding = true(num_mov_nodes,1);
    
    % Loop through each time step creating a correspondence between moving
    % and reference surface values
    for i = 1 : t_steps
        % Get the current configuration of the model
        x = NodeArray + u{i};
        xs_mov = x(ns_mov,:);
        
        % Get the closest point projection between surfaces
        %if i == 1
        cppStruct = ClosestPointTriSurfV2( fa_ref, x, xs_mov(is_sliding,:) );
        %end
        dist = cppStruct.dist;
        inside = cppStruct.inside;
        
        % Identify sliding points
        is_sliding_new = (dist < dist_tol) | inside;
        f_cp = cppStruct.face(is_sliding_new);
        Ncp = cppStruct.Ncp(is_sliding_new,:);
        xcp = cppStruct.Xcp(is_sliding_new,:);
        is_sliding(is_sliding) = is_sliding_new; 
        num_pairs = sum(is_sliding);
        np_mov = ns_mov(is_sliding);
        xp_mov = xs_mov(is_sliding,:);
        
        for j = 1:num_pairs         
            % Get the corresponding reference location for the cp
            Xcp{i}(np_mov(j),:) = Ncp(j,:) * NodeArray(fa_ref(f_cp(j),:),:);
            
            % Get the corresponding displacements for the cp
            Ucp{i}(np_mov(j),:) = Ncp(j,:) * u{end}(fa_ref(f_cp(j),:),:);
            
            if i < t_steps
                % Get the incremental displacement for the moving point
                u_mov_inc = u{i+1}(np_mov(j),:) - u{i}(np_mov(j),:);

                % Get the corresponding incremental displacement for the cp
                ucp_inc = Ncp(j,:) * u{i+1}(fa_ref(f_cp(j),:),:)...
                          - Ncp(j,:) * u{i}(fa_ref(f_cp(j),:),:);

                % Get the relative displacement between the two points
                u_rel = u_mov_inc - ucp_inc;
                
                % Define the normal vector as the difference vector
                % between the moving point and the projected point
                normal = xp_mov(j,:) - xcp(j,:);
                normal = normal / norm(normal);

                % Get the tangential component of the relative displacement
                % between the two points
                ut_inc{i}(np_mov(j),:) = u_rel - (normal*u_rel')*normal;
                %ut_inc{i}(np_mov(j),:) = u_rel;
            end            
        end
    end
    
%% Clean up and derive final results
    % Initialize final arrays
    traj = cell(1,t_steps);
    traj{1} = NodeArray(np_mov,:);
    s_total = zeros(size(np_mov,1),1);
    s_total_conv = zeros(size(np_mov,1),1);
    
    for i = 1 : t_steps
        % Remove non-sliding nodes from Xcp and Ucp
        Xcp{i} = Xcp{i}(np_mov,:);
        Ucp{i} = Ucp{i}(np_mov,:);
        
        % Store spatial trajectory and cumulative sliding distances
        if i > 1
            traj{i} = traj{i-1} + ut_inc{i-1}(np_mov,:); 
            s_total = s_total + vecnorm( ut_inc{i-1}(np_mov,:), 2, 2 );
            s_total_conv = s_total_conv + vecnorm( Xcp{i} - Xcp{i-1}, 2, 2 );
        end
    end
    
    traj_alt = cell2mat(traj);
    Xcp_alt = cell2mat(Xcp);

%% Plot Spatial Sliding Trajectories Against Moving Surface Mesh
%     % To reduce figure clutter sample to a given spacing
%     Xp_mov = NodeArray(np_mov,:);
%     if plot_blank > 0
%         n_samp = GridSample(Xp_mov,plot_blank);
%     else
%         n_samp = 1 : size(Xp_mov,1);
%     end
%     
%     % Plot sliding trajectories against mesh
%     if ~plot_hide
%         figure()
%         hold on
%         for i = 1:length(n_samp)
%             a = traj_alt(n_samp(i),:);
%             a = reshape(a,3,[]);
%             a = a';
%             % Apply a scale factor to trajectory
%             a = a*plot_scale - a(1,:)*(plot_scale - 1);
%             plot3(a(:,1),a(:,2),a(:,3),'LineWidth',2)
%         end
%         plot3(Xp_mov(n_samp,1),Xp_mov(n_samp,2),Xp_mov(n_samp,3),'.k','MarkerSize',15)
%         trimesh(fa_mov,NodeArray(:,1),NodeArray(:,2),NodeArray(:,3), ...
%             'EdgeColor',[0.5,0.5,0.5])
%         if exist('subject','var')
%             title(sprintf('Spatial Sliding Trajectories: %s \nBlack dot indicates starting point of trajectory',subject))
%         else
%             title(sprintf('Spatial Sliding Trajectories \nBlack dot indicates starting point of trajectory'))
%         end
% 
%         % Camera control
%         daspect([1 1 1])
%         if plot_invert
%             set(gca, 'Zdir', 'reverse')
%             set(gca, 'Ydir', 'reverse')
%         end    
%         hold off
%     end

%% Plot convected sliding trajectory and contour plots together
    % Get s_total for the full NodeArray
    s_total_full = nan( size(NodeArray,1), 1 );
    s_total_full(np_mov) = s_total;
    
    % To reduce figure clutter sample to a given spacing
    Xcp_start = Xcp{1};
    if plot_blank > 0
        n_samp = GridSample(Xcp_start,plot_blank);
    else
        n_samp = 1 : size(Xcp_start,1);
    end
    
    if ~plot_hide
        figure()        
        % Plot convected sliding trajectories against mesh
        sp_array(1) = subplot(1,2,1);
        hold on
        for i = 1:length(n_samp)
            a = Xcp_alt(n_samp(i),:);
            a = reshape(a,3,[]);
            a = a';
            plot3(a(:,1),a(:,2),a(:,3),'LineWidth',2)
        end
        plot3(Xcp_start(n_samp,1),Xcp_start(n_samp,2),Xcp_start(n_samp,3),'.k','MarkerSize',15)
        trimesh(fa_ref,NodeArray(:,1),NodeArray(:,2),NodeArray(:,3), ...
            'EdgeColor',[0.5,0.5,0.5])
        % Camera control
        daspect([1 1 1])
        if plot_invert
            set(gca, 'Zdir', 'reverse')
            set(gca, 'Ydir', 'reverse')
        end
        title(sprintf('Convected Sliding Trajectories \nBlack dot indicates starting point of trajectory'))
        hold off
        
        % Plot sliding distance contour
        sp_array(2) = subplot(1,2,2);
        patch('Faces',fa_mov,'Vertices',NodeArray,'FaceVertexCData',s_total_full,...
            'FaceColor','interp','EdgeColor','None')
        % Colorbar control
        colorbar()
        caxis([ min(s_total), max(s_total) ])
        %caxis([ 0, 30 ])
        % Camera Control
        daspect([1 1 1])
        if plot_invert
            set(gca, 'Zdir', 'reverse')
            set(gca, 'Ydir', 'reverse')
        end
        title(sprintf('Spatial Sliding Distance Contour'))
        
        % Link camera between subplots
        hlink = linkprop( sp_array, {'CameraPosition','CameraUpVector'} );
    end
    
%% Generate outStruct
    outStruct.NodeArray = NodeArray;
    outStruct.ElementArray = ElementArray;
    outStruct.nodes_mov = np_mov;
    outStruct.U_mov = u{end}(np_mov,:);
    outStruct.Xcp = Xcp;
    outStruct.Ucp = Ucp;
    outStruct.s_total = s_total;
    outStruct.s_total_conv = s_total_conv;
    outStruct.traj = traj;
    if ~plot_hide
        outStruct.hlink = hlink;
    end
end