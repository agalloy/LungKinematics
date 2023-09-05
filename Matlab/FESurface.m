

function SurfaceArray = FESurface(ElementArray)
    %Put all of the outward faces of the mesh into an array
    FA = [ElementArray(:,[3 2 1]);ElementArray(:,[2 3 4]);ElementArray(:,[3 1 4]);ElementArray(:,[4 1 2])];
    %Sort the FaceArray so order of nodes doesn't matter
    FA_sorted = sort(FA,2);
    %Get each unique face
    [FA_unique,sort_index,unique_index] = unique(FA_sorted,'rows');    
    %Count the number of times each face appears in mesh
    FaceCount = histc(unique_index,1:length(FA_unique));
    %Find the faces that only appear once (surface faces)
    SurfaceArray = FA(sort_index(FaceCount == 1),:);
end