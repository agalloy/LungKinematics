% Plots vector field streamlines on a triangulated surface

function PlotTriSurfStreamline( FaceArray, NodeArray, X, options )
%% Parse options structure
    num_dim = size(NodeArray,2);
    num_face = size(FaceArray,1);
    % Get mean edge length for defaults
    edge_vec = [ NodeArray(FaceArray(:,2),:) - NodeArray(FaceArray(:,1),:)
                 NodeArray(FaceArray(:,3),:) - NodeArray(FaceArray(:,2),:)
                 NodeArray(FaceArray(:,1),:) - NodeArray(FaceArray(:,3),:) ];
    edge_mag = vecnorm(edge_vec,2,2);
    mean_length = mean(edge_mag);
             
    if isfield(options,'seeds')
        seeds = options.seeds;
    elseif isfield(options,'spacing')
        spacing = options.spacing;
        if num_dim == 3
            seeds = GridSample( NodeArray, spacing );
        else
            error('Spacing option is only supported for 3D surfaces')
        end
    else
        if num_dim == 3
            spacing = 6;
            seeds = GridSample( NodeArray, spacing );
        end
        if num_dim == 2
            seeds = (1:10:size(NodeArray,1))';
        end
    end
    if isfield(options,'num_steps')
        num_steps = options.num_steps;
    else
        num_steps = ceil( sqrt(num_face) );
    end
    if isfield(options,'step_size')
        step_size = options.step_size;
    else
        step_size = mean_length / 2;
    end
    if isfield(options,'LineWidth')
        LineWidth = options.LineWidth;
    else
        LineWidth = 2;
    end
    if isfield(options,'min_mag')
        min_mag = options.min_mag;
    else
        min_mag = step_size * 10^-10;
    end
        
%% Trace streamlines for each seed
    num_seeds = size(seeds,1);
    % Initialize important arrays
    % Cell array for storing streamlines
    stream = cell( num_seeds, 1 );
    % Sx3 array for storing the position of each stream at the current time
    % step
    Pos = NodeArray(seeds,:);
    % Sx3 array for storing the interpolated value of X at the current time
    % step
    Xi = X(seeds,:);
    % Algorithm for 3D surfaces
    if num_dim == 3
        for i = 1:num_steps
            % Get vector magnitude
            X_mag = vecnorm(Xi,2,2);            
            % Ignore small vectors to avoid numerical issues
            ignore = X_mag < min_mag;
            % Take a step in the direction of the vector field
            Pos3D = Pos;
            Pos3D(~ignore,:) = Pos3D(~ignore,:) + Xi(~ignore,:)./X_mag(~ignore)*step_size;
            % Reproject the current position back onto the surface
            cppStruct = ClosestPointTriSurfV2( FaceArray, NodeArray, Pos3D );
            Pos = cppStruct.Xcp;

            % Interpolate the value of X at the new position
            face = cppStruct.face;
            Ncp = cppStruct.Ncp; % N x I 
            for j = 1:num_seeds
                % Get X at each face vertex
                Xf = X(FaceArray(face(j),:),:); % N x X
                Xi(j,:) = Ncp(j,:) * Xf;

                % Update the new positions in stream
                stream{j}(i,:) = Pos(j,:);
            end
        end
    % Algorithm for 2D surfaces
    elseif num_dim == 2
        % Create a triangulation object of the 2D surface
        TR = triangulation(FaceArray,NodeArray);
        for i=1:num_steps
            % Get vector magnitude
            X_mag = vecnorm(Xi,2,2);     
            % Ignore small vectors to avoid numerical issues
            ignore = X_mag < min_mag;
            % Take a step in the direction of the vector field
            Pos_new = Pos;
            Pos_new(~ignore,:) = Pos(~ignore,:) + Xi(~ignore,:)./X_mag(~ignore)*step_size;
            
            for j = 1:num_seeds
                % Get the traingle associated with the current location
                % barycentric coordinates at that location
                face = pointLocation( TR, Pos_new(j,:) );
                if ~isnan(face)
                    Pos(j,:) = Pos_new(j,:);
                    Ncp = cartesianToBarycentric( TR, face, Pos(j,:) );
                    % Get the value of X at each vertex of that triangle
                    Xf = X(FaceArray(face,:),:);
                    % Interpolate to the current seed position
                    Xi(j,:) = Ncp * Xf;
                end               
                % Update the new positions in stream
                stream{j}(i,:) = Pos(j,:);
            end
        end
    end
%% Plot streamline trajectories
    hold on
    if num_dim == 3
        trimesh( FaceArray, NodeArray(:,1), NodeArray(:,2), NodeArray(:,3), ...
            'FaceColor', [0.80,0.80,0.80], 'EdgeColor', [0,0,0], 'FaceAlpha', 1 )
        for i = 1:num_seeds
            plot3( stream{i}(:,1),stream{i}(:,2),stream{i}(:,3), ...
                'Color', [0,0,1], 'LineWidth', LineWidth )
        end
        daspect([1 1 1])
    elseif num_dim == 2
        trimesh( FaceArray, NodeArray(:,1), NodeArray(:,2), zeros(size(NodeArray,1),1), ...
            'FaceColor', [0.80,0.80,0.80], 'EdgeColor', [0,0,0], 'FaceAlpha', 1 )
        for i = 1:num_seeds
            plot3( stream{i}(:,1),stream{i}(:,2), zeros(num_steps,1), ...
                'Color', [0,0,1], 'LineWidth', LineWidth )
        end
        daspect([1 1 1])
    end
end