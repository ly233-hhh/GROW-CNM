clc
clear
tic
% Set to your local example data path
workdir = fullfile(pwd, 'example_data'); 

% Example weight map (e.g., a T-statistic map)
weight_map_path = fullfile(workdir, 'Twosample2_t1_tstat2.nii'); 
mni_coord = [-19.5, -6, -12];

% Target cluster size in voxels. (e.g., 206459 mm^3 divided by voxel volume 1.5^3)
target_size = round(206459 / (1.5^3)); 

output_path = fullfile(workdir, 'generated_cluster.nii');

tic
generate_cluster(weight_map_path, mni_coord, target_size, output_path);
toc