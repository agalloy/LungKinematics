% Display the original vector field, its three HHD components, and the
% (co)potential.
% Input:
%   NodeArray = N x 3 or N x 2 array of node positions
%   FaceArray = E x 3 array with face connectivity info
%   HHD_struct = The output structure of the HHD_GradientRecon function
%   options = (options) Additional options

%% Main Function
function hlink = ViewHHDComponents( FaceArray, NodeArray, HHD_struct, options )
    % Parse optional options input
    if ~exist( 'options', 'var' )
        options = struct();
    end
    % Plot invert option
    if isfield(options,'plot_invert')
        plot_invert = options.plot_invert;
    else
        plot_invert = false;
    end     

% Load important data
    n_dim = size(NodeArray,2);
    title_string = ["Potential","Original Vector Field","Copotential"
                    "Exact Component","Harmonic Component","Coexact Component"]';
    plot_data = {HHD_struct.alpha, HHD_struct.omega(:,1:n_dim), HHD_struct.beta_n
                 HHD_struct.diff_alpha(:,1:n_dim), HHD_struct.gamma(:,1:n_dim), HHD_struct.codiff_beta(:,1:n_dim) }';
    scalar = str2func('PlotScalarField');
    vector = str2func('PlotVectorField');
    plot_type = {scalar, vector, scalar
                 vector, vector, vector }';

% Loop through each of the six plots to make
    figure()
    for i = 1:6
        sp_array(i) = subplot(2,3,i);
        hold on
        title(title_string(i))
        % Create the correct type of plot for the data
        plot_type{i}( FaceArray, NodeArray, plot_data{i}, options );
        % Set correct aspect ratio
        daspect([1 1 1])        
        % Flip plot if desired
        if plot_invert
            set(gca, 'Zdir', 'reverse')
            set(gca, 'Ydir', 'reverse')
        end
        hold off
    end
    % Link the camera of the plots
    hlink = linkprop( sp_array, {'CameraPosition','CameraUpVector'} );
end

%% Function to plot a scalar contour plot
function PlotScalarField( FaceArray, NodeArray, A, ~ )
    patch( 'Faces', FaceArray, 'Vertices', NodeArray, 'FaceColor', 'interp', 'CData', A, ...
           'EdgeAlpha', 0 );
    colorbar()
end

%% Function to plot a vector streamline plot
function PlotVectorField( FaceArray, NodeArray, X, options )
    PlotTriSurfStreamline( FaceArray, NodeArray, X, options )
end