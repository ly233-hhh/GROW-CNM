function generate_cluster(weight_map_path, mni_coord, target_size, output_path)
% GENERATE_CLUSTER Generates a fixed-size irregular cluster based on MNI coordinates and a weight map (SPM).
% 
% Algorithm: Global Optimal Region Growing (Best-First Region Growing)
% =========================================================================
% This algorithm maintains a global candidate queue of all unassigned boundary voxels.
% In each iteration, it adds the voxel with the highest weight from the global queue to the cluster.
% This ensures the cluster grows along the highest statistical values, perfectly aligning with 
% the true high-value clusters in the weight map.
% =========================================================================
%
% Usage:
%   generate_cluster(weight_map_path, mni_coord, target_size, output_path)
%
% Inputs:
%   weight_map_path : Path to the weight map NIfTI file (e.g., t-statistic map).
%   mni_coord       : Seed MNI coordinates [x, y, z].
%   target_size     : Target cluster size in number of voxels.
%   output_path     : Output binary NIfTI mask path.

    if nargin < 4
        error('Insufficient inputs. Usage: generate_cluster(weight_map_path, mni_coord, target_size, output_path)');
    end
    if numel(mni_coord) ~= 3
        error('mni_coord must be a vector of length 3, e.g., [24 -10 36].');
    end
    if target_size <= 0 || round(target_size) ~= target_size
        error('target_size must be a positive integer.');
    end

    if exist('spm_vol', 'file') ~= 2 || exist('spm_read_vols', 'file') ~= 2 || exist('spm_write_vol', 'file') ~= 2
        error('SPM functions not found. Please add SPM12 to your MATLAB path.');
    end

    fprintf('1. Reading NIfTI file using SPM...\n');
    V = spm_vol(weight_map_path);
    weight_map = spm_read_vols(V);

    if numel(V) > 1
        warning('4D NIfTI detected. Only the first volume will be used.');
        V = V(1);
        weight_map = weight_map(:, :, :, 1);
    end

    dim = size(weight_map);
    if numel(dim) > 3
        dim = dim(1:3);
        weight_map = weight_map(:, :, :, 1);
    end

    invalid_weight = ~isfinite(weight_map);
    if any(invalid_weight(:))
        weight_map(invalid_weight) = -Inf;
    end

    fprintf('2. Converting MNI coordinates to voxel coordinates...\n');
    mni_col = [mni_coord(:); 1];
    voxel_col = inv(V.mat) * mni_col;
    seed_voxel = round(voxel_col(1:3))';

    if any(seed_voxel < 1) || any(seed_voxel > dim)
        error('Seed voxel coordinate out of bounds. Check MNI or affine matrix.');
    end

    fprintf('3. Initializing global optimal queue, target size: %d voxels...\n', target_size);
    cluster_mask = false(dim);
    visited = false(dim);

    seed_idx = sub2ind(dim, seed_voxel(1), seed_voxel(2), seed_voxel(3));
    cluster_mask(seed_idx) = true;
    visited(seed_idx) = true;

    [dx, dy, dz] = ndgrid(-1:1, -1:1, -1:1);
    neighbors_offset = [dx(:), dy(:), dz(:)];
    neighbors_offset(~any(neighbors_offset, 2), :) = [];

    boundary_idx = [];
    current_size = 1;
    [boundary_idx, visited] = add_neighbors_spm(seed_voxel, dim, visited, neighbors_offset, boundary_idx, weight_map);

    fprintf('4. Starting best-first region growing...\n');
    while current_size < target_size && ~isempty(boundary_idx)
        b_weights = weight_map(boundary_idx);
        [best_score, best_idx_in_queue] = max(b_weights);
        
        if isinf(best_score) && best_score < 0
            warning('All remaining candidates have invalid weights, stopping early.');
            break;
        end

        best_voxel_idx = boundary_idx(best_idx_in_queue);
        cluster_mask(best_voxel_idx) = true;
        current_size = current_size + 1;
        boundary_idx(best_idx_in_queue) = [];

        [vx, vy, vz] = ind2sub(dim, best_voxel_idx);
        [boundary_idx, visited] = add_neighbors_spm([vx, vy, vz], dim, visited, neighbors_offset, boundary_idx, weight_map);

        if mod(current_size, 1000) == 0 || current_size == target_size
            fprintf('   Generated %d / %d voxels...\n', current_size, target_size);
        end
    end

    fprintf('5. Saving results...\n');
    Vout = V;
    Vout.fname = output_path;
    Vout.dt = [spm_type('uint8') 0];
    Vout.pinfo = [1; 0; 0];
    spm_write_vol(Vout, uint8(cluster_mask));
    fprintf('Done.\n');
end

function [boundary_idx, visited] = add_neighbors_spm(voxel, dim, visited, neighbors_offset, boundary_idx, weight_map)
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