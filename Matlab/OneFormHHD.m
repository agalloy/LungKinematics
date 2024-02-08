% A function that performs Helmholtz-Hodge Decomposition directly on
% 1-forms (instead of taking a vector field and converting it to one).
% Inputs:
%   FaceArray = F X 3 node connectivity array
%   NodeArray = N X 3 array of node positions
%   omega = E X 1 a discrete 1-form
%   options = A struct of additional options
% Outputs:
%   alpha = N X 1 scalar (0-form) potential
%   beta = F X 1 2-form copotential

function [alpha,beta] = OneFormHHD(FaceArray,NodeArray,omega,options)
%% Parse options structure or apply defaults
    % Parse optional options input
    if ~exist( 'options', 'var' )
        options = struct();
    end
    % Boundary conditions to use
        % bc = 0: Potential is set to zero at a single boundary node 
        %   (Harmonic component is absorbed into exact component)
        % bc = 1: Potential is set to zero at all boundary nodes 
        %   (Harmonic component is orthogonal to exact and coexact components)
    if isfield(options,'bc')
        bc = options.bc;
    else
        bc = 1;
    end     
    % Display verification info
    if isfield(options,'verify')
        verify = options.verify;
    else
        verify = false;
    end     
    % Pre-computed geometry info
    if isfield(options,'DEC')
        DEC = options.DEC;
    else
        DEC = AssembleDEC(FaceArray,NodeArray);
    end
    
    % Load geometric prerequisites
    b_nodes = DEC.b_nodes;
    d0 = DEC.d0;
    d1 = DEC.d1;
    hs1 = DEC.hs1;
    hs2 = DEC.hs2;
    
    num_nodes = size(NodeArray,1);
    num_faces = size(FaceArray,1);

%% Compute the potential function (0-form) and its associated 1-form
    % Set boundary conditions
    alpha_bc = nan(num_nodes,1);
    if bitand(bc,1) == 1
        bc_nodes = b_nodes;
        % Set potential on all boundary nodes to 0
        alpha_bc(bc_nodes) = 0;
    else
        % Set potential to zero at a single point
        bc_nodes = (1:num_nodes)' == find(b_nodes,1);
        alpha_bc(bc_nodes) = 0;
    end
    
    % Compute "stiffness matrix"
    K = d0' * hs1 * d0;
    
    % Compute right hand side vector (with BC's)
    f = d0' * hs1 * omega;
    f_bc = K(~bc_nodes,bc_nodes) * alpha_bc(bc_nodes);
    RHS = f(~bc_nodes) - f_bc;
    
    % Solve system for alpha
    alpha = nan(num_nodes,1);
    alpha(bc_nodes) = alpha_bc(bc_nodes);
    alpha(~bc_nodes) = K(~bc_nodes,~bc_nodes) \ RHS;  

%% Compute the copotential function (2-form) and its associated 1-form
    % Set boundary conditions
    beta_bc = nan(num_faces,1);
    bc_faces = false(num_faces,1);
    % No boundary conditions defined for copotential as of now
    % However, discretization automatically enforces 0 copotential along boundary 
    
    % Compute "stiffness matrix"
    K = d1 * hs1^(-1) * d1' * hs2;
    
    % Compute right hand side vector (with BC's)
    f = d1 * omega;
    f_bc = K(~bc_faces,bc_faces) * beta_bc(bc_faces);
    RHS = f(~bc_faces) - f_bc;
    
    % Solve system for beta
    beta = nan(num_faces,1);
    beta(bc_faces) = beta_bc(bc_faces);
    beta(~bc_faces) = K(~bc_faces,~bc_faces) \ RHS;
    
%% Verification info
    if verify
        % Residual for each solved quantity
        alpha_res = d0' * hs1 * omega - d0' * hs1 * d0 * alpha;
        disp('Max residual: Potential')
        disp( max( alpha_res(~bc_nodes) ) )
        beta_res = d1 * omega - d1 * hs1^(-1) * d1' * hs2 * beta;
        disp('Max residual: Copotential')
        disp( max( beta_res(~bc_faces) ) )
    end
end