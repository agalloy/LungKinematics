

function vel = GenerateVelocityField( NodeArray, options )
%% Parse options
    if isfield(options,'t_type')
        t_type = options.t_type;
    else
        t_type = 'identity';
    end
    
    if strcmp( t_type, 'identity' )
        arg1 = [];
        arg2 = [];
    elseif strcmp( t_type, 'translation' )
        arg1 = options.t_vel;
        arg2 = [];
    elseif strcmp( t_type, 'rotation' )
        arg1 = options.ang_vel;
        arg2 = options.origin;
    elseif strcmp( t_type, 'stretch' )
        arg1 = options.scale_factor;
        arg2 = options.origin;
    else
        error('This transformation type is not supported!')
    end
    
%% Call the correct function based on inputs
    % Evaluate the vector field at each node
    handle = str2func( t_type );
    vel = handle(NodeArray,arg1,arg2);
    
end

%% Local Functions
% Identity VF
function V = identity(X,~,~) 
    V = X.*0;
end

% Translation VF
function V = translation(X,k,~)
    V = X.*0 + k;
end

% Rotation VF
function V = rotation(X,k,o)
    V = cross( X - o, repmat(k,size(X,1),1) );
end

% Stretch VF
function V = stretch(X,k,o)
    V = (X-o) .* (k-[1,1,1]);
end