function coordout=coord_mni2vox(coordin,opt,T)
% COORD_MNI2VOX Coordinate transform between MNI152 and voxel space
% Input:
%   coordin : Nx3 Matrix
%   opt     : 'm2v' = mni->voxel, 'v2m' = voxel->mni
%   T       : affine matrix, optional
% Output:
%   coordout: output coordinates, Nx3

if nargin<1
    usage;
    return
end
t0= [ -2     0     0    92;...
      0     2     0  -128;...
      0     0     2   -74;...
      0     0     0     1];

if nargin<2
    opt='m2v';
end
if nargin <3
    T=t0;
end
if isempty(T)
    T=t0;
end
if isempty(opt)
    opt='m2v';
end

if strcmp(opt,'m2v')
    coordout = [coordin(:,1) coordin(:,2) coordin(:,3) ones(size(coordin,1),1)]*(inv(T))';
    coordout(:,4) = [];
    coordout = round(coordout);
elseif strcmp(opt,'v2m')
    coordout = T*[coordin(:,1) coordin(:,2) coordin(:,3) ones(size(coordin,1),1)]';
    coordout = coordout';
    coordout(:,4) = [];
end
end

function usage()
fprintf('====Coordinate transform between MNI152 and voxel space====\n');
fprintf('Format: coordout=coord_mni2vox(coordin,opt,T)\n');
fprintf('Input:\n--coordin: input coordinates, Nx3 Matrix\n');
fprintf('--opt: ''m2v''=mni->voxel, ''v2m''=voxel->mni\n');
fprintf('--T: affine matrix, optional\n');
fprintf('Output:\n--coordout: output coordinates, Nx3\n');
end