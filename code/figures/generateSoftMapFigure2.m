clear;
files = {'SHREC_281_285_straight','SHREC_281_285_expdiff','SHREC_281_285_reverse'};%'SHREC_281_285_expdiff'

sources = [ 293 ... % snout
            246 ... % left bicep
            906 ... % right knee
            387 ... % midsection
            941 ... % left thigh
            361  ... % left ankle
            535 ... % right forearm
            ];

for i=1:length(files)
    load(sprintf('../data/maps/%s.mat',files{i}));
    close all
    sphereRadius = .05;
    f = map(:,sources);
    generateSoftMapMeshes(M1,M2,'sourcemap.obj',[files{i} '.obj'],fixed1,fixed2,sources,f,sphereRadius);
end