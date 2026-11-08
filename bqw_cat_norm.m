function ydefs=bqw_cat_norm(data,options)
if nargin<1
    fprintf('===========normalize using cat12============\n');
    fprintf('Usage:ydefs=bqw_cat_norm(data,options)\n--data: image used to segment & normalize\n');
    fprintf('--options.voxsize: output voxel size\n--options.interp: interpolation method,0=nn,1=trilinear,2-4=spline\n');
    fprintf('--options.prefix: output prefix\no--options.bbox: output bonding box\n');
    fprintf('----------Written by QinWen 20250618-------------\n');
    return
end
cat=cat_get_defaults;
tpmpath=cat.extopts.darteltpm{1,1};
if nargin<2
    options.bbox=spm_get_bbox(tpmpath);
    vol=spm_vol(tpmpath);
    options.voxsize=abs(eig(vol(1).mat(1:3,1:3)))';
    options.interp=1;
    options.prefix='mni_';
end
if ~isfield(options,'bbox')
    options.bbox=spm_get_bbox(tpmpath);
end
if ~isfield(options,'voxsize')
    vol=spm_vol(tpmpath);
    options.voxsize=abs(eig(vol(1).mat(1:3,1:3)))';
end
if ~isfield(options,'interp')
    options.interp=1;
end
if ~isfield(options,'prefix')
    options.prefix='mni_';
end


[indir,fname,ext]=fileparts(data);
if isempty(indir)
    indir=pwd;
end
indir=[indir,filesep,'mri'];
ydefs{1,1}=[indir,filesep,'y_',fname,'.nii'];
ydefs{2,1}=[indir,filesep,'iy_',fname,'.nii'];
if ~exist(ydefs{1,1},'file')
    cat_job(data,options);
    fprintf('..Done cat segmentation: %s\n',data);
end
outmT1=[indir,filesep,'..',filesep,options.prefix,'m',fname,ext];
if ~exist(outmT1,'file')
    copyfile([indir,filesep,'wm',fname,ext],outmT1);
    fprintf('..rename output T1: %s\n',outmT1);
end
% warp other images
if isfield(options,'otherimages')
    for f=1:length(options.otherimages)
        other=options.otherimages{f};
        [odir,oname,oext]=fileparts(other);
        outother=[odir,filesep,options.prefix,oname,oext];
        if ~exist(outother,'file')
            vols{1}=[other,',1'];
            dartel_write(vols,ydefs,options);
            fprintf('..done normalize: %s\n',other);
        end
    end
end
end
%% subfunction: dartel_write
function dartel_write(data,ydefs,options)
matlabbatch{1}.spm.spatial.normalise.write.subj.def = ydefs(1);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = data;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = options.bbox;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = options.voxsize;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = options.interp;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = options.prefix;
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
end

function cat_job(data,opt)
if nargin<2
    cat.extopts.vox=1.5;
else
    cat.extopts.vox=opt.voxsize(1);
end
%%default parameters
cat.extopts.WMHC=0; % segment WMHC
cat.extopts.SLC=0; % segment stroke lesion
cat.output.surface=0;
cat.output.ROImenu.noROI = struct([]);
cat.output.sROImenu.noROI = struct([]);
cat.output.ROI = 0;
cat.extopts.nproc=1;
% bias T1
cat.output.bias.native=1;
cat.output.bias.warped=1;
cat.output.bias.dartel=0;
% label
cat.output.label.native=0;
cat.output.label.warped=0;
cat.output.bias.dartel=0;
% Las (local bais corrected)
cat.output.las.native=1;
cat.output.las.warped=0;
cat.output.las.dartel=0;
%GM
cat.output.GM.native=1;
cat.output.GM.warped=0;
cat.output.GM.mod=0;
cat.output.GM.dartel=0;
% WM
cat.output.WM.native=1;
cat.output.WM.warped=0;
cat.output.WM.mod=0;
cat.output.WM.dartel=0;
% CSF
cat.output.CSF.native=1;
cat.output.CSF.warped=0;
cat.output.CSF.mod=0;
cat.output.CSF.dartel=0;
%WMH
cat.output.WMH.native=0;
cat.output.WMH.warped=0;
cat.output.WMH.mod=0;
cat.output.WMH.dartel=0;
%stroke lesion
cat.output.SL.native=0;
cat.output.SL.warped=0;
cat.output.SL.mod=0;
cat.output.SL.dartel=0;
% TPMC for dartel
cat.output.TPMC.native=0;
cat.output.TPMC.warped=0;
cat.output.TPMC.mod=0;
cat.output.TPMC.dartel=0;
% y_ deformations [forward backward]
cat.output.warps = [1 1];
% jacobian determinant 0/1 (none/yes)
cat.output.jacobian.warped = 1;

%% main
matlabbatch{1}.spm.tools.cat.estwrite = cat;
matlabbatch{1}.spm.tools.cat.estwrite.data = {data};
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
end
