% Identifies the contact interface between two surfaces
% Inputs:
%   F_mov = Fm X 3 Face connectivity array for moving surface
%   N_mov = Nm X 3 Node coordinate array (may contain nodes not on the
%       surface as long as numbering is consistent with F_mov) 
%   F_ref = Fr X 3 Face connectivity array for reference surface
%   N_ref = Nr X 3 Node coordinate array (may contain nodes not on the
%       surface as long as numbering is consistent with F_ref) 
%   options = (Optional) Structure with additional options
% Outputs:
%   faces_ci = An array of face identifiers for the moving surface
%       indicating faces on the contact interface
%   nodes_ci = An array of node identifiers for the moving surface
%       indicating nodes on the contact interface

function [ FID_ci, NID_ci ] = ContactInterface( FA_mov, NA_mov, FA_ref, NA_ref, options )
%% Parse options structure
    if ~exist( 'options', 'var' )
        options = struct();
    end
    % Read options or apply defaults
    if isfield(options,'dist_tol')
        dist_tol = options.dist_tol;
    else
        dist_tol = 0.5;
    end    

%% Get the closest point projection between surfaces
    % Query only the nodes in NA_mov that appear in FA_mov
    NID_mov = unique( FA_mov );
    projStruct = ClosestPointTriSurfV2( FA_ref, NA_ref, NA_mov(NID_mov,:) );
        
    % Determine the sliding points
    dist = projStruct.dist;
    inside = projStruct.inside;
    in_contact = (dist < dist_tol) | (inside & (dist < 10*dist_tol));
    
%% Return outputs
    NID_ci = NID_mov( in_contact );
    FID_ci = find(all( ismember(FA_mov,NID_ci), 2 ));

end