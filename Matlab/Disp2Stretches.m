% Get principal stretches from displacement data. Works for linear
% tetrahedral elements. Also return element volumes for conveniece with 
% averaging.

function [lambda,volume] = Disp2Stretches(NodeArray,ElementArray,DispArray)
%% Loop through each element computing the principal stretches
    % Initialize Arrays
    lambda = nan( size(ElementArray,1), 3 );
    volume = nan( size(ElementArray,1), 1 );
    for i = 1:size(ElementArray,1)
        X = NodeArray( ElementArray(i,:), : );
        u = DispArray( ElementArray(i,:), : );
        
        % Get the deformation gradient
        M =     [X(1,1), X(2,1), X(3,1), X(4,1)
                 X(1,2), X(2,2), X(3,2), X(4,2)
                 X(1,3), X(2,3), X(3,3), X(4,3)
                 1,      1,      1,      1];
        dN_dX = M^(-1);
        dN_dX = dN_dX(:,1:3);
        du_dX = u' * dN_dX;
        F = eye(3) + du_dX;
        
        % Get eigenvalues from right Cauchy Green to compute principal stretches
        C = F'*F;
        eig_C = sort( eig(C), 'descend' );
        lambda(i,:) = sqrt(eig_C');
        
        % Get the volume of the tetrahedron
        volume(i) = abs(det(M))/6;
    end

end