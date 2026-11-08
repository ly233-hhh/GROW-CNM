function [coordvox, coordvec, mask] = SphereMask(center, radius,vol,opt)
% initialize
if nargin<3
    usage
    return
end
if nargin<4
    opt.center_type='mm';
    opt.omask='';
    opt.brainmask=[];
    opt.display='no';
end
if ~isfield(opt,'display')
    opt.display='no';
end
if ~isfield(opt,'center_type')
    opt.center_type='mm';
end
if ~isfield(opt,'omask')
    opt.omask='';
end
if ~isfield(opt,'brainmask')
    opt.brainmask=[];
end
matrixSize=vol.dim;
T=vol.mat;
% convert mm coords to voxel unit
if strcmp(opt.center_type,'mm')
    center=coord_mm2vox(center,'m2v',T);
end
% convert mm radius to voxel unit
voxelSize = sqrt(sum(T(1:3,1:3).^2, 1)); % 
radius_vox=radius/mean(voxelSize); % divided by mean voxel size

% generate sphere mask
mx = matrixSize(1);
my = matrixSize(2);
mz = matrixSize(3);
[X, Y, Z] = ndgrid(1:mx, 1:my, 1:mz);

% generate sphere ball for each center and merge
ncoord=size(center,1);
mask=zeros(vol.dim);
for n =1:ncoord
    x0 = center(n,1);
    y0 = center(n,2);
    z0 = center(n,3);
    dist = sqrt((X - x0).^2 + (Y - y0).^2 + (Z - z0).^2);
    % generate sphere mask
    mask0 = dist <= radius_vox;
    mask0=double(mask0);
    mask = mask+mask0;
end
mask(mask>0)=1; %binary overlapped sphere voxels
% remove voxels outside the brain
if ~isempty(opt.brainmask)
    if isnumeric(opt.brainmask)
        brainmask=opt.brainmask;
    elseif exist(opt.brainmask,'file')
        volm=spm_vol(opt.brainmask);
        imgm=spm_read_vols(volm);
        brainmask=imgm;
    else
        % consider refimage as brain mask
        volm=vol;
        imgm=spm_read_vols(volm);
        brainmask=imgm;
        brainmask(brainmask~=0)=1;
    end
    masksize=size(brainmask);
    if ~any(masksize==matrixSize)
        error('un-matched dimension between brainmask and vol.dim\n');
    end
    mask=mask.*brainmask;
end

% generate voxel indices vector
coordvec = find(mask);
[x, y, z] = ind2sub([mx, my, mz], coordvec);
coordvox = [x, y, z];

% display the balls
if strcmp(opt.display,'yes')
    figure;
    p0 = patch(isosurface(brainmask, 0.5));
    isonormals(brainmask, p0);
    set(p0, 'FaceColor',[0.1 0.1 0.1], 'EdgeColor', 'none', 'FaceAlpha', 0.05);
    hold on
    p = patch(isosurface(mask, 0.5));
    isonormals(mask, p);
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    daspect([1 1 1]);
    camlight; lighting gouraud;
    view(-90,90);
    axis off
end

% save mask image and mask figure
if ~isempty(opt.omask)
    [outdir,outname]=fileparts(opt.omask);
    if isempty(outdir)
        outdir=pwd;
    end
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    vol.fname=opt.omask;
    vol.dt=[2 0];
    spm_write_vol(vol,mask);
    if strcmp(opt.display,'yes')
        outfig=[outdir,filesep,outname,'.fig'];
        saveas(gcf,outfig);
    end
    hold off
end
fprintf('-------------Done sphere mask generatation ( Written by QinWen 20250621)------------\n');
end

%% subfunctions
function coordout=coord_mm2vox(coordin,opt,T)
% coordin: N*3 coordinates
% opt: 'm2v'=mm to voxel conversion,'v2m' voxel to mm conversion
% T: affine matrix
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

% usage
function usage()
fprintf('--------------Generate sphere mask--------------\nInput:\n');
fprintf('--center: sphere center coordinates [ncoord,3],mm or voxel coordinates determined by opt.center_type \n')
fprintf('--radius: sphere radius(mm)\n');
fprintf('--vol: spm_vol which define the mask space\n');
fprintf('--opt.center_type:center coordinate type,''mm'',''vox''(default)\n');
fprintf('--opt.brainmask: mask image to remove non-brain voxels,file path or matrix\n')
fprintf('--opt.display: display the sphere ROI,''yes'',''no''(default)\n')
fprintf('--opt.omask: mask image output path,''''=do not save mask image\n');
fprintf('Output:\n--coordvox: voxel coordinates [Nvox, 3]\n--coordvec: voxel index vector\n--mask: mask matrix\n]');
fprintf('-----------Written by Qin Wen at 20250621------------\n');
end