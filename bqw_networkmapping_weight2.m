function bqw_networkmapping_weight2(seed, connectome_folders, OutSubjDir, options)
%% 初始化参数
t0 = clock;
if nargin < 3
    usage;
    return
end
if nargin < 4
    options = struct();
end

% 基本默认值
def_options.conn_type = 'FN';
def_options.conn_pref = 'fcGSRdwB1aroma_fmap';
def_options.conn_smooth = 6;
def_options.conn_npca = 10;
def_options.seed_type = 'coord';
def_options.seed_coord_weight = 'virtualcluster';
def_options.seed_weight_scale = 'none';
def_options.seed_image_weight = 'all';
def_options.seed_coord_radius = 6;
def_options.brainmask = '';
def_options.seed_restrict_mask = '';   % 已废弃，不再使用
def_options.auto_switch_to_image = false;
def_options.auto_switch_image_weight = 'all';
def_options.feature_collinear_thr = 0.95;
def_options.rm_individual = 'yes';
def_options.seed_coord_oname = 'coord_virtual_seed';
def_options.seed_cluster_unit = 'mm';
def_options.weight_map = '';           % 新增：虚拟簇生长依赖的权重图路径

fnames = fieldnames(def_options);
for i = 1:length(fnames)
    if ~isfield(options, fnames{i})
        options.(fnames{i}) = def_options.(fnames{i});
    end
end

[~, subjname] = fileparts(OutSubjDir);
[indir, ~] = fileparts(connectome_folders{1});

fprintf('--------Calculate Network Mapping: %s --------\n', subjname);
fprintf('..Connectome Folder: %s\n', indir);
fprintf('..Connectome Type: %s \n', options.conn_type);
fprintf('..Connectome File Pref: %s\n', options.conn_pref);
fprintf('..Seed Type: %s\n', options.seed_type);
fprintf('..Seed Coord Weight: %s\n', options.seed_coord_weight);
fprintf('..Seed weight scale: %s\n', options.seed_weight_scale);
fprintf('..Brain Mask: %s\n', options.brainmask);
if ~isempty(options.seed_restrict_mask)
    fprintf('..Restrict Mask: %s (ignored for virtualcluster modes)\n', options.seed_restrict_mask);
end
fprintf('..Output Subject Directory: %s\n', OutSubjDir);

if ~exist(OutSubjDir, 'dir')
    mkdir(OutSubjDir);
end

%% Step 1: 确定参考图像和脑模板
fprintf('Step 1: Define reference image and mask\n');
if strcmp(options.conn_type, 'SN')
    refimage = [connectome_folders{1}, filesep, options.conn_pref, '_density.nii'];
    if exist(refimage, 'file')
        volref = spm_vol(refimage);
    else
        error('--Missing fiber density image: %s\n', refimage);
    end
else
    tmp = dir([connectome_folders{1}, filesep, options.conn_pref, '*.nii*']);
    if isempty(tmp)
        error('--Missing images in folder(%s) with prefix(%s)\n', connectome_folders{1}, options.conn_pref);
    else
        refimage = [connectome_folders{1}, filesep, tmp(1).name];
        volref = spm_vol([refimage, ',1']);
    end
end
fprintf('..Reference image is: %s\n', refimage);

if ~isempty(options.brainmask)
    if exist(options.brainmask, 'file')
        volmask = spm_vol(options.brainmask);
        [maskdir, maskname, maskext] = fileparts(options.brainmask);
        if any(volmask.mat(:) ~= volref.mat(:)) || any(volmask.dim ~= volref.dim)
            args.interp = 0;
            args.pref = 'resl_';
            args.overwrite = 0;
            rmask = fullfile(maskdir, [args.pref, maskname, maskext]);
            ormask = fullfile(OutSubjDir, [args.pref, maskname, maskext]);
            if ~exist(ormask, 'file')
                bqw_reslice(options.brainmask, refimage, args);
                movefile(rmask, ormask);
            end
            options.brainmask = ormask;
            volmask = spm_vol(options.brainmask);
            fprintf('..Resliced Brain Mask is: %s\n', options.brainmask);
        end
    else
        warning('..Mask image missing, omit: %s', options.brainmask);
    end
end

%% Step 2: 种子预处理（混合模式）
fprintf('Step 2: Seed Preprocessing \n');
if ~strcmp(options.seed_type, 'coord')
    error('This function only supports seed_type=''coord''.');
end

fprintf('..Seed is a study''s coordinates group \n');
[ncoord, ndim] = size(seed);
if ndim < 3
    error('seed coordinates should have at least 3 column [x y z stats clustersize]');
end
if ndim < 4, seed(:,4) = 1; end
if ndim < 5, seed(:,5) = 1; end

weight_mode = options.seed_coord_weight;
if strcmp(weight_mode, 'all')
    modes_to_calc = {'equal', 'stats', 'virtualcluster', 'both'};
else
    modes_to_calc = {weight_mode};
end

need_sphere = any(ismember(modes_to_calc, {'equal', 'stats'}));
need_virtual = any(ismember(modes_to_calc, {'virtualcluster', 'both'}));

if need_virtual
    if any(seed(:,5) <= 1)
        error('Virtual cluster requires cluster_size > 1 for all seeds.');
    end
    if isempty(options.brainmask) || ~exist(options.brainmask, 'file')
        error('Virtual cluster requires a valid brain mask file.');
    end
    % 不再要求 restrict_mask
    if isempty(options.weight_map) || ~exist(options.weight_map, 'file')
        error('Virtual cluster (region growing) requires a weight_map file (e.g., t-map).');
    end
end

% ---------- 球形种子 (equal/stats) ----------
if need_sphere
    fprintf('....Generating sphere-based seeds (for equal/stats)\n');
    coordimgout = fullfile(OutSubjDir, [options.seed_coord_oname, '_sphere_binary.nii']);
    statsimgout = fullfile(OutSubjDir, [options.seed_coord_oname, '_sphere_stats.nii']);
    
    % 若任一文件缺失，则重新生成（避免残存文件导致报错）
    if ~exist(coordimgout, 'file') || ~exist(statsimgout, 'file')
        volseed = volref;
        opt.center_type = 'mm';
        opt.brainmask = options.brainmask;
        opt.display = 'no';
        opt.omask = '';
        centers = seed(:,1:3);
        seed_coordmat_sphere = zeros(volseed.dim);
        seed_stats_map = zeros(volseed.dim);
        for n = 1:ncoord
            [cvox, cvec, ~] = SphereMask(centers(n,:), options.seed_coord_radius, volseed, opt);
            
            % ★ 修复：若stats值为NaN或空，设为默认值1，确保球体能被标记
            current_stat = seed(n,4);
            if isnan(current_stat) || isempty(current_stat)
                current_stat = 1;   % 或其他正数，保证球体存在
            end
            
            for i = 1:length(cvec)
                idx = cvec(i);
                if abs(current_stat) > abs(seed_stats_map(idx))
                    seed_stats_map(idx) = current_stat;
                    seed_coordmat_sphere(idx) = 1;
                end
            end
        end
        seed_coordvec_sphere = find(seed_coordmat_sphere);
        [x,y,z] = ind2sub(volseed.dim, seed_coordvec_sphere);
        seed_coordvox_sphere = [x,y,z];
        volimg = volseed; volimg.dt = [16,0];
        volimg.fname = coordimgout;
        spm_write_vol(volimg, seed_coordmat_sphere);
        volimg.fname = fullfile(OutSubjDir, [options.seed_coord_oname, '_sphere_stats.nii']);
        spm_write_vol(volimg, seed_stats_map);
    else
        volseed = spm_vol(coordimgout);
        seed_coordmat_sphere = spm_read_vols(volseed);
        seed_coordvec_sphere = find(seed_coordmat_sphere);
        [x,y,z] = ind2sub(volseed.dim, seed_coordvec_sphere);
        seed_coordvox_sphere = [x,y,z];
        vol1 = spm_vol(fullfile(OutSubjDir, [options.seed_coord_oname, '_sphere_stats.nii']));
        seed_stats_map = spm_read_vols(vol1);
    end
    stats_sphere = seed_stats_map(seed_coordvec_sphere);
    switch options.seed_weight_scale
        case 'sqrt'
            stats_sphere = sign(stats_sphere) .* sqrt(abs(stats_sphere));
        case 'log'
            stats_sphere = sign(stats_sphere) .* log(abs(stats_sphere) + 1);
    end
    stats_abs_sphere = abs(stats_sphere);

    sphere_weights = [];
    sphere_weight_names = {};
    if ismember('equal', modes_to_calc)
        sphere_weights = [sphere_weights, ones(length(seed_coordvec_sphere), 1)];
        sphere_weight_names{end+1} = 'wEqual';
    end
    if ismember('stats', modes_to_calc)
        sphere_weights = [sphere_weights, stats_abs_sphere];
        sphere_weight_names{end+1} = ['w', options.seed_weight_scale, 'Stats'];
    end
else
    sphere_weights = [];
    sphere_weight_names = {};
    seed_coordvec_sphere = [];
    seed_coordvox_sphere = [];
end

% ---------- 虚拟簇种子 (virtualcluster/both) 基于权重图的区域生长 ----------
if need_virtual
    fprintf('....Generating virtual cluster seeds via weight-map-based region growing\n');
    volseed = volref;

    % ---------- 读取权重图（若空间不匹配则自动重采样）----------
    Vw = spm_vol(options.weight_map);
    if any(Vw.mat(:) ~= volseed.mat(:)) || any(Vw.dim ~= volseed.dim)
        fprintf('....Weight map space mismatch, reslicing to reference image space...\n');
        args.interp = 1;              % 1 = 线性插值，适合连续统计图
        args.pref   = 'resl_';        % 输出文件前缀
        args.overwrite = 0;           % 若已存在则不重复生成
        [wdir, wname, wext] = fileparts(options.weight_map);
        rweight = fullfile(wdir, [args.pref, wname, wext]);  % 重采样后的完整路径
        if ~exist(rweight, 'file')
            bqw_reslice(options.weight_map, refimage, args);
        end
        options.weight_map = rweight;
        Vw = spm_vol(rweight);
        fprintf('....Resliced weight map: %s\n', rweight);
    end
    weight_map_data = spm_read_vols(Vw);
    weight_map_data(~isfinite(weight_map_data)) = -Inf;

    % 体素体积，用于 mm³ -> 体素数转换
    voxel_volume = abs(det(volseed.mat));

    % 初始化掩模和权重图
    seed_coordmat_v = zeros(volseed.dim);
    weight_img = zeros(volseed.dim);

    % ★ 关键修改：维护一个全局已占据掩膜，保证团块不重叠
    occupied_mask = false(volseed.dim);

    % 对每个种子点独立生长，但需避开前面种子已占据的区域
    for n = 1:ncoord
        mni_seed = seed(n, 1:3);
        cluster_size = seed(n, 5);   % 单位由 seed_cluster_unit 决定

        if strcmp(options.seed_cluster_unit, 'mm')
            target_vox = round(cluster_size / voxel_volume);
        else % 'vox'
            target_vox = round(cluster_size);
        end
        if target_vox < 1
            warning('Seed %d cluster size too small after conversion, skipping.', n);
            continue;
        end

        fprintf('....Growing cluster %d/%d, MNI=[%.1f,%.1f,%.1f], target vox=%d\n', ...
                n, ncoord, mni_seed, target_vox);

        % 将已占据的体素从权重图中排除，防止重叠
        temp_weight = weight_map_data;
        temp_weight(occupied_mask) = -Inf;

        % 调用生长函数，使用临时权重图
        cluster_mask_i = grow_cluster_from_weight(temp_weight, mni_seed, target_vox, volseed);

        % 更新占据掩膜
        occupied_mask = occupied_mask | cluster_mask_i;

        % 合并到总掩膜
        seed_coordmat_v = seed_coordmat_v | cluster_mask_i;

        % 若需要 both 加权，保存权重图上的原始值（不重叠部分）
        if ismember('both', modes_to_calc)
            % 只对新生成的团块区域赋值，且保留原始权重（未缩放）
            weight_img(cluster_mask_i) = weight_map_data(cluster_mask_i);
        end
    end

    % 若只要 virtualcluster（等权重），则将权重全部设为1
    if ~ismember('both', modes_to_calc)
        weight_img(seed_coordmat_v > 0) = 1;
    end

    % 不再应用 restrict_mask，直接保存
    volimg = volseed; volimg.dt = [16,0];
    volimg.fname = fullfile(OutSubjDir, [options.seed_coord_oname, '_virtualcluster_binary.nii']);
    spm_write_vol(volimg, seed_coordmat_v);
    volimg.fname = fullfile(OutSubjDir, [options.seed_coord_oname, '_virtualcluster_weight.nii']);
    spm_write_vol(volimg, weight_img);

    seed_coordvec_virtual = find(seed_coordmat_v);
    if isempty(seed_coordvec_virtual)
        error('No virtual seed voxels generated; check weight map or cluster size.');
    end
    [x,y,z] = ind2sub(size(seed_coordmat_v), seed_coordvec_virtual);
    seed_coordvox_virtual = [x,y,z];

    stats_virtual = weight_img(seed_coordvec_virtual);
    switch options.seed_weight_scale
        case 'sqrt'
            stats_virtual = sign(stats_virtual) .* sqrt(abs(stats_virtual));
        case 'log'
            stats_virtual = sign(stats_virtual) .* log(abs(stats_virtual) + 1);
    end
    stats_abs_virtual = abs(stats_virtual);

    % 构建虚拟簇组的权重
    virtual_weights = [];
    virtual_weight_names = {};
    if ismember('virtualcluster', modes_to_calc)
        virtual_weights = [virtual_weights, ones(length(seed_coordvec_virtual), 1)];
        virtual_weight_names{end+1} = 'wVirtualCluster';
    end
    if ismember('both', modes_to_calc)
        virtual_weights = [virtual_weights, stats_abs_virtual];
        virtual_weight_names{end+1} = ['w', options.seed_weight_scale, 'StatsVirtualCluster'];
    end
else
    virtual_weights = [];
    virtual_weight_names = {};
    seed_coordvec_virtual = [];
    seed_coordvox_virtual = [];
end

%% Step 3: 网络映射计算（分组处理，完全保留原有逻辑）
switch options.conn_type
    case 'FN'
        fprintf('Step 3: Calculate Individual FN-based CNM (grouped modes)\n');
        groups = {};
        if need_sphere && ~isempty(sphere_weights)
            groups{end+1} = struct('name', 'Sphere', ...
                'coordvec', seed_coordvec_sphere, ...
                'coordvox', seed_coordvox_sphere, ...
                'weights', sphere_weights, ...
                'weight_names', {sphere_weight_names});
        end
        if need_virtual && ~isempty(virtual_weights)
            groups{end+1} = struct('name', 'Virtual', ...
                'coordvec', seed_coordvec_virtual, ...
                'coordvox', seed_coordvox_virtual, ...
                'weights', virtual_weights, ...
                'weight_names', {virtual_weight_names});
        end

        if isempty(groups)
            error('No valid seed group to calculate.');
        end

        if exist('volmask','var')
            mask = spm_read_vols(volmask);
            mask2D = mask(:) > 0;
        else
            mask2D = [];
        end

        for g = 1:length(groups)
            thisGroup = groups{g};
            fprintf('   --- Processing group: %s (weights: %s) ---\n', ...
                thisGroup.name, strjoin(thisGroup.weight_names, ', '));
            seed_coordvec_i = thisGroup.coordvec;
            seed_weight_mat = thisGroup.weights;
            weight_names_i = thisGroup.weight_names;
            nWeights = size(seed_weight_mat, 2);

            FCdir = fullfile(OutSubjDir, 'FC');
            if ~exist(FCdir, 'dir'), mkdir(FCdir); end

            all_exist = true;
            for ww = 1:nWeights
                mFCpath = fullfile(OutSubjDir, ['mean_zFC_Map_', weight_names_i{ww}, '.nii']);
                if ~exist(mFCpath, 'file')
                    all_exist = false;
                    break;
                end
            end
            if all_exist
                fprintf('   ..All mean/T maps already exist for this group, skipping.\n');
                continue;
            end

            for s = 1:length(connectome_folders)
                tic
                [~, connsub] = fileparts(connectome_folders{s});
                tmp = dir(fullfile(connectome_folders{s}, [options.conn_pref, '*.nii*']));
                if isempty(tmp)
                    warning('--Missing fMRI images in %s\n', connectome_folders{s});
                    continue;
                end
                volimg = [];
                for f = 1:length(tmp)
                    tmppath = fullfile(connectome_folders{s}, tmp(f).name);
                    if f == 1
                        volimg = spm_vol(tmppath);
                    else
                        volimg = [volimg; spm_vol(tmppath)];
                    end
                end
                fmri = spm_read_vols(volimg);
                fmri(isnan(fmri)) = 0;
                fmri(isinf(fmri)) = 0;
                msize = size(fmri);
                fmri2D = reshape(fmri, prod(msize(1:3)), msize(4));
                if ~isempty(mask2D)
                    fmri2D(~mask2D, :) = 0;
                end
                vid = (max(fmri2D, [], 2) - min(fmri2D, [], 2)) > 0;
                vox_tp = fmri2D(vid, :)';

                tp_mat = fmri2D(seed_coordvec_i, :)' * seed_weight_mat;
                for ww = 1:nWeights
                    outpref = fullfile(FCdir, [connsub, '_zFC']);
                    FC_file = [outpref, '_', weight_names_i{ww}, '.nii'];
                    if exist(FC_file, 'file')
                        fprintf('....FC already exist: %s\n', FC_file);
                        continue;
                    end
                    tp_w = tp_mat(:, ww);
                    r = corr(vox_tp, tp_w);
                    zFC = 0.5 * log((1 + r) ./ (1 - r));
                    zFCout = zeros(msize(1:3));
                    zFCout(vid) = zFC;
                    volout = volimg(1);
                    volout.dt = [16, 0];
                    volout.fname = FC_file;
                    spm_write_vol(volout, zFCout);
                    fprintf('....Done FC for %s, %s\n', connsub, weight_names_i{ww});
                end
                toc
            end

            % 组水平统计
            for ww = 1:nWeights
                wname = weight_names_i{ww};
                FClist = dir(fullfile(FCdir, ['*_zFC_', wname, '.nii']));
                if isempty(FClist)
                    warning('No individual FC maps for %s', wname);
                    continue;
                end
                szFCfile = fullfile(OutSubjDir, ['s', num2str(options.conn_smooth), '_zFC_', wname, '.nii']);
                if ~exist(szFCfile, 'file')
                    for f = 1:length(FClist)
                        volfc = spm_vol(fullfile(FCdir, FClist(f).name));
                        outsm = volfc;
                        outsm.fname = szFCfile;
                        outsm.n(1) = f;
                        spm_smooth(volfc, outsm, options.conn_smooth);
                    end
                end
                vols = spm_vol(szFCfile);
                imgfc = spm_read_vols(vols);
                msize = size(imgfc);
                if exist('volmask','var')
                    mask = spm_read_vols(volmask);
                else
                    mask = spm_read_vols(vols(1)) ~= 0;
                end
                mid = mask(:) ~= 0;
                imgfc2d = reshape(imgfc, prod(msize(1:3)), msize(4));
                vid = std(imgfc2d, [], 2) > 0;
                fcid = vid & mid;
                Y0 = imgfc2d(fcid, :)';
                dfe = msize(4) - 1;
                [mzFC, stdFC] = grpstats(Y0, ones(dfe+1,1), {'mean', 'std'});
                Tval = mzFC * sqrt(msize(4)) ./ stdFC;
                mean_zFC = zeros(msize(1:3));
                mean_zFC(fcid) = mzFC;
                Tval_zFC = zeros(msize(1:3));
                Tval_zFC(fcid) = Tval';
                volout = vols(1);
                volout.descrip = sprintf('dfe=%d;smooth=%dmm', dfe, options.conn_smooth);
                volout.fname = fullfile(OutSubjDir, ['mean_zFC_Map_', wname, '.nii']);
                spm_write_vol(volout, mean_zFC);
                volout.fname = fullfile(OutSubjDir, ['Tvals_zFC_Map_', wname, '.nii']);
                spm_write_vol(volout, Tval_zFC);
                if strcmp(options.rm_individual, 'yes')
                    delete(fullfile(FCdir, ['*_zFC_', wname, '.nii']));
                    delete(szFCfile);
                end
                fprintf('....Done Tval/mean zFC Maps for %s\n', wname);
            end
        end

    otherwise
        error('Only FN mode is fully adapted. Extend SN/TSN similarly.');
end

t1 = clock;
et = etime(t1, t0);
fprintf('-------- Finish Network Mapping: %s ---------\n Total Spend %0.1f seconds.\n\n\n\n', subjname, et);
end

%% ==============================================
%% 新增子函数：基于权重图的全局最优区域生长
function cluster_mask = grow_cluster_from_weight(weight_map, mni_coord, target_size, V_ref)
% 从MNI坐标出发，在权重图上使用优先队列（Best-First）生长至指定体素数
    dim = V_ref.dim;
    mni_col = [mni_coord(:); 1];
    voxel_col = V_ref.mat \ mni_col;
    seed_voxel = round(voxel_col(1:3)');
    if any(seed_voxel < 1) || any(seed_voxel > dim)
        error('Seed coordinate outside image bounds.');
    end

    cluster_mask = false(dim);
    visited = false(dim);

    seed_idx = sub2ind(dim, seed_voxel(1), seed_voxel(2), seed_voxel(3));
    cluster_mask(seed_idx) = true;
    visited(seed_idx) = true;
    current_size = 1;

    % 26邻域偏移
    [dx, dy, dz] = ndgrid(-1:1, -1:1, -1:1);
    neighbors_offset = [dx(:), dy(:), dz(:)];
    neighbors_offset(~any(neighbors_offset, 2), :) = [];

    boundary_idx = [];
    [boundary_idx, visited] = add_boundary_neighbors(seed_voxel, dim, visited, neighbors_offset, boundary_idx, weight_map);

    while current_size < target_size && ~isempty(boundary_idx)
        b_weights = weight_map(boundary_idx);
        [~, best_idx_in_queue] = max(b_weights);
        if isinf(b_weights(best_idx_in_queue))
            warning('Reached invalid weight region, stopping early.');
            break;
        end
        best_voxel_idx = boundary_idx(best_idx_in_queue);
        cluster_mask(best_voxel_idx) = true;
        current_size = current_size + 1;
        boundary_idx(best_idx_in_queue) = [];
        [vx, vy, vz] = ind2sub(dim, best_voxel_idx);
        [boundary_idx, visited] = add_boundary_neighbors([vx, vy, vz], dim, visited, neighbors_offset, boundary_idx, weight_map);
        if mod(current_size, 10000) == 0
            fprintf('      ...%d/%d voxels\n', current_size, target_size);
        end
    end
end

function [boundary_idx, visited] = add_boundary_neighbors(voxel, dim, visited, neighbors_offset, boundary_idx, weight_map)
    nx = voxel(1) + neighbors_offset(:, 1);
    ny = voxel(2) + neighbors_offset(:, 2);
    nz = voxel(3) + neighbors_offset(:, 3);
    valid_mask = nx >= 1 & nx <= dim(1) & ...
                 ny >= 1 & ny <= dim(2) & ...
                 nz >= 1 & nz <= dim(3);
    nx = nx(valid_mask); ny = ny(valid_mask); nz = nz(valid_mask);
    if isempty(nx), return; end
    n_idx = sub2ind(dim, nx, ny, nz);
    unvisited_mask = ~visited(n_idx) & isfinite(weight_map(n_idx));
    new_idx = n_idx(unvisited_mask);
    boundary_idx = [boundary_idx; new_idx(:)];
    visited(new_idx) = true;
end

%% ==============================================
%% 其余子函数（原样保留，未修改）
function seed_fiber_density(fiberin,refimage,seed_coordvox,seed_weight,seed_weight_name,outpref,options)
% ...（保持原函数完整）...
if nargin<6
    fprintf('----Calculate fiber density (nfiber/mm3) passing through the seed----\n');
    fprintf('Format: seed_fiber_density(fiberin,refimage,seed_coordvox,seed_weight,seed_weight_name,outpref,options)\n');
    fprintf('Input:\n--fiberin: input fibers, ftr/dsistudio .mat\n');
    fprintf('--refimage: image space define the density map\n');
    fprintf('--seed_coordvox: seed coordinates in voxel\n--seed_weight: seed voxel weights\n');
    fprintf('--seed_weight_name: seed weight name {nweight,1}\n--outpref: output density image directory and pref\n');
    return
end
tic
[fibdir,fibname,fibext]=fileparts(fiberin);
if strcmp(fibext,'.trk')
    fib_coords=bmio_read_trk(fiberin);
    fibers=vertcat(fib_coords.coord);
    idx=vertcat(fib_coords.nPoints);
    fibers(:,4) = repelem(1:length(idx), idx)';
else
    load(fiberin);
end
vol=spm_vol(refimage);
msize=vol.dim;
fibdensity_img=zeros(msize);
voxelsize=sqrt(sum(vol.mat(1:3,1:3).^2, 1));
radius=max(voxelsize);
seedfiber= SeedFiberRound(fibers,seed_coordvox,radius,options.brainmask);
if  strcmp(options.seed_image_weight,'equal') && strcmp(options.seed_type,'image')
    seed_fib_weight=ones(size(seedfiber,1),1);
elseif strcmp(options.seed_coord_weight,'equal') && strcmp(options.seed_type,'coord')
    seed_fib_weight=ones(size(seedfiber,1),1);
else
    seed_fib_weight=define_fib_weights(seedfiber,seed_coordvox,seed_weight);
end
fibcoordvx=ceil(seedfiber);
fibcoordvx(fibcoordvx(:,1)>msize(1)|fibcoordvx(:,1)<1,:)=[];
fibcoordvx(fibcoordvx(:,2)>msize(2)|fibcoordvx(:,2)<1,:)=[];
fibcoordvx(fibcoordvx(:,3)>msize(3)|fibcoordvx(:,3)<1,:)=[];
[coords,ia]=unique(fibcoordvx,'rows');
seed_fib_weight_uniq=seed_fib_weight(ia,:);
[coorduniq,~,ic] = unique(coords(:,1:3),'rows');
for w=1:size(seed_weight,2)
    output=[outpref,'_',seed_weight_name{w},'.nii'];
    wfd = accumarray(ic,seed_fib_weight_uniq(:,w));
    fibdensity_coords = [coorduniq, wfd];
    for i=1:size(fibdensity_coords,1)
        xyz=fibdensity_coords(i,1:3);
        fibdensity_img(xyz(1),xyz(2),xyz(3))=fibdensity_coords(i,4);
    end
    voxvol=abs(det(vol.mat));
    fibdensity_img=fibdensity_img/voxvol;
    vol.fname=output;
    vol.dt=[16 0];
    spm_write_vol(vol,fibdensity_img);
end
toc
end

function seedfiber = SeedFiberRound(fibers,seed_coordvox,radius,brainmask)
points=round(fibers(:,1:3)/radius)*radius;
if exist(brainmask,'file')
    volm=spm_vol(brainmask);
    mask=spm_read_vols(volm);
    msize=size(mask);
    maskvox=coord2scalar(find(mask(:)),msize,'s2c');
    mid=ismember(points,maskvox,'rows');
else
    mid=ones(size(points,1),1);
end
centers=seed_coordvox;
sid=ismember(points,centers,'rows');
fibidx=fibers(:,4);
uniq_seedfibid=unique(fibidx(sid));
fid=ismember(fibidx,uniq_seedfibid);
seedfiber=fibers(fid&mid,:);
fprintf('......Done Extract fibers pass the seed: \n');
end

function seed_fib_weight = define_fib_weights(seedfiber,seed_coordvox,seed_weight)
ns = KDTreeSearcher(seed_coordvox);
[idx_seed, dist] = knnsearch(ns, seedfiber (:, 1:3));
data = [seedfiber(:, 4), dist, idx_seed];
nfibcoord=accumarray(seedfiber(:,4),1);
nfibcoord(nfibcoord==0)=[];
fibidx=unique(seedfiber(:, 4));
seed_fib_weight=zeros(size(seedfiber,1),size(seed_weight,2));
for n = 1:length(fibidx)
    i0=sum(nfibcoord(1:n-1))+1;
    i1=sum(nfibcoord(1:n));
    fibid=i0:i1;
    tmpdata=data(fibid,:);
    [~,minidx]=min(tmpdata(:,2));
    seed_idx=tmpdata(minidx,3);
    seed_fib_weight(fibid,:)=repmat(seed_weight(seed_idx,:),length(fibid),1);
end
fprintf('......Done define fiber weights: \n');
end

function usage ()
fprintf('==============Individualized Network Mapping===============\n')
fprintf('--seed:seed for a subject, can be a nifti mask file for LNM/ANM, or [n 3] coordinates for CNM\n')
fprintf('--connectome_folders: normative connectome folders[rootdir,filesep,subid], cell string(one line per individual)\n');
fprintf('--OutSubjDir:output directory for seed subject\n');
fprintf('--options.brainmask: REQUIRED brain mask file to remove fake connectome and voxels; defines reference space\n');
fprintf('--options.seed_cluster_unit: cluster size unit,''mm''(default) or ''vox''\n');
fprintf('--options.conn_type: connectome type, ''SN''=structural network (SDC),''FN''=functional network (FDC),''TSN''=texture similarity network \n');
fprintf('--options.conn_pref: prefix for connectome file\n');
fprintf('--options.conn_smooth: smooth kernel of individual connectome map (t-test need in FC/TSN), 6(mm)\n');
fprintf('--options.seed_type: seed type,''image''=seed mask image,''coord''=coordinates,[x,y,z,stats,clustersize]\n')
fprintf('--options.seed_image_weight: seed voxel weights,''equal''=equal weight,''abs''=weight by absolute values\n');
fprintf('                          ''all''=output all the above .\n');
fprintf('--options.seed_coord_weight: seed coord weights,''equal''=equal weight,''cluster''=weight by cluster size,''virtualcluster''=generate virtualcluster by coordinate cluster size\n');
fprintf('                          ''stats''=weight by statistics (T/F),''both''=weight by cluster size and statistics,''all''=output all the above \n');
fprintf('--options.seed_weight_scale: weight scaling, ''none''(default),''sqrt'',''log''')
fprintf('--options.seed_coord_radius: sphere radius (mm) of coords,6(default)\n')
fprintf('--options.seed_coord_oname: filename of output coordiniates mask\n')
fprintf('--options.rm_individual: delete individual norm connectome ,''yes'' (default)\n')
fprintf('--options.weight_map: (NEW) path to weight map for virtual cluster growing.\n')
fprintf('--------------Written by QinWen at 20250619----------------\n')
fprintf('--------------Revised by QinWen at 20251024, added virtualcluster generation with the same size as the seed coordniates----------------\n')
fprintf('--------------Revised 2026, replaced VirtualClusterMask with weight-map-based region growing. Removed restrict_mask requirement.----------------\n')
end
