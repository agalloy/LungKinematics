% Compute the surface area of a triangulated surface.
% Inputs:
%   FaceArray = [F x 3] Face connectivity array
%   NodeArray = [N x 3] Node position array
% Outputs:
%   A_total = Surface area of entire surface
%   A_tri = Surface area of individual triangles

function [A_total, A_tri] = TriSurfArea(FaceArray, NodeArray)
    e12 = NodeArray(FaceArray(:,2),:) - NodeArray(FaceArray(:,1),:);
    e13 = NodeArray(FaceArray(:,3),:) - NodeArray(FaceArray(:,1),:);
    A_tri = vecnorm( cross(e12,e13,2), 2, 2 )/2;
    A_total = sum(A_tri);
end