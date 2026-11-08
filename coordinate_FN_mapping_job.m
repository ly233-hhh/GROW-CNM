function coordinate_FN_mapping_job
% COORDINATE_FN_MAPPING_JOB 
% Example script to run GROW-CNM mapping for a given seed file.

clc; clear;

%% -------- Read environment variables or set paths --------
seedfile = getenv('SEEDFILE');
if isempty(seedfile) || exist(seedfile, 'file') ~= 2
    % 本地测试：可直接在这里修改路径
    seedfile = 'path/to/your/seed_data.txt'; 
end

outdir_root = getenv('OUTROOT');
if isempty(outdir_root)
    outdir_root = 'path/to/output_directory';
end
[~, study] = fileparts(seedfile);

% Weight map (e.g., T-statistic map used for GROW region growing)
weight_map_path = 'path/to/weight_map.nii'; % 【请修改为真实路径】

% Connectome subject list and root path
conn_sublist    = 'path/to/subject_list.txt'; % 包含受试者ID的文件
connectome_root = 'path/to/connectome_root/'; 

% Brain mask
brainmask_path = 'path/to/brainmask.nii';

%% -------- Check files --------
if ~exist(seedfile, 'file'), error('Seed file not found: %s', seedfile); end
if ~exist(weight_map_path, 'file'), error('Weight map not found: %s', weight_map_path); end
if ~exist(brainmask_path, 'file'), error('Brain mask not found: %s', brainmask_path); end

OutStudyDir = fullfile(outdir_root, study);

%% -------- Set options --------
options.conn_type          = 'FN';
options.conn_pref          = 'fc_fmap'; 
options.conn_smooth        = 6;               
options.seed_type          = 'coord';         
options.seed_coord_weight  = 'virtualcluster'; % 使用 GROW 算法
options.seed_weight_scale  = 'none';          
options.seed_cluster_unit  = 'mm';            
options.seed_coord_radius  = 6;               
options.rm_individual      = 'yes';           
options.brainmask          = brainmask_path;
options.weight_map         = weight_map_path; 

%% -------- Add paths --------
addpath('path/to/core_functions'); % 核心函数库目录
% addpath('path/to/spm12'); spm('defaults', 'fmri');

%% -------- Read seeds --------
fprintf('Reading seed file: %s\n', seedfile);
if endsWith(seedfile, {'.xlsx', '.xls'})
    seedcoord = xlsread(seedfile);
else
    seedcoord = load(seedfile);
end
if size(seedcoord, 2) < 5
    error('Seed coordinates need at least 5 columns: [x y z stats clustersize]');
end

%% -------- Build connectome folders --------
fid = fopen(conn_sublist, 'r');
conn_subid = textscan(fid, '%s'); fclose(fid);
connectome_folders = fullfile(connectome_root, conn_subid{1});

%% -------- Run core analysis --------
fprintf('Starting network mapping...\n');
bqw_networkmapping_weight2(seedcoord, connectome_folders, OutStudyDir, options);
fprintf('Done. Results saved in: %s\n', OutStudyDir);