clear;
files = {'SCAPE_00_05_straight','SCAPE_00_05_expdiff','SCAPE_00_05_reverse'};

sources = [ 984 ... % face
            736 ... % left bicep
            233 ... % right knee
            511 ... % midsection
            329 ... % left thigh
            92  ... % left ankle
            653 ... % right forearm
            ];

for i=1:length(files)
    load(sprintf('../data/maps/%s.mat',files{i}));
    close all
    sphereRadius = .05;
    f = map(:,sources);
    generateSoftMapMeshes(M1,M2,'sourcemap.obj',[files{i} '.obj'],fixed1,fixed2,sources,f,sphereRadius);
end