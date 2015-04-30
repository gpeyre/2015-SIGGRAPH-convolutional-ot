clear;
files = {'SHREC_388_392_straight'};

sources = [ 196 ... % back right leg top
            527 ... % mid back
            691 ... % front right leg
            87  ... % front left knee
            980 ... % forehead
            529 ... % mid-side
            971 ... %ear
            ];

for i=1:length(files)
    load(sprintf('../data/maps/%s.mat',files{i}));
    close all
    sphereRadius = .05;
    f = map(:,sources);
    generateSoftMapMeshes(M1,M2,'sourcemap3.obj',[files{i} '.obj'],fixed1,fixed2,sources,f,sphereRadius);
end