function E2 = subdivideGraphEdges(E)

E2 = [];

ne = max(E(:));
for i = 1:size(E,1)
    v1 = E(i,1);
    v2 = E(i,2);
    vNew = ne+i;
    E2 = [E2 ; v1 vNew ; vNew v2];
end