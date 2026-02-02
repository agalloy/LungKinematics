

function [V_total, V_elem] = TetVolume( ElementArray, NodeArray )
    e12 = NodeArray(ElementArray(:,2),:) - NodeArray(ElementArray(:,1),:);
    e13 = NodeArray(ElementArray(:,3),:) - NodeArray(ElementArray(:,1),:);
    e14 = NodeArray(ElementArray(:,4),:) - NodeArray(ElementArray(:,1),:);
    V_elem = abs(dot( cross(e12,e13,2), e14, 2 )/6);
    V_total = sum(V_elem);
end