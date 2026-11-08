function bqw_reslice(images,reslvox,opt)
if nargin<1
    fprintf('\n========= Reslice Nifti Images =========\n');
    fprintf('Format: bqw_reslice(images,reslvox,opt)\n');
    fprintf('--images: source images need resliced, textlist or cellarray\n');
    fprintf('--reslvox: resliced voxel, can be voxelsize (numeric), or reference image path(string)\n');
    fprintf('--opt.interp: interporation method,0=''neighbor'',1=''trilinear'',2-4=''spline''\n');
    fprintf('--opt.pref: prefix of resliced images\n--opt.bb:bonding box\n');
    fprintf('--opt.overwrite: 1=''yes'',0=''no'',\n----------Written by QinWen at 20201019-------------\n');
    return
end

if nargin<2
    reslvox=[1 1 1];
end

if nargin<3
    opt.interp=1;
    opt.pref='resl_';
    opt.overwrite=0;
    opt.bb=[];
end

if ~isfield(opt,'interp')
    opt.interp=1;
end
if ~isfield(opt,'pref')
    opt.pref='resl_';
end
if ~isfield(opt,'overwrite')
    opt.overwrite=0;
end
if ~isfield(opt,'bb') % add bondingbox optional inputs20201203
    opt.bb=[];
end

if ischar(images)
    if exist(images,'file')
        [tmpdir,fname,ext]=fileparts(images);
        if strcmp(ext,'.gz') || strcmp(ext,'.nii')%correct at 20210225
            imgpaths={images};
        else
            fid=fopen(images);
            tmp=textscan(fid,'%s');
            imgpaths=tmp{1};
            fclose(fid);
        end
    else
        error('missing file path %s\n',images);
    end
elseif iscellstr(images)
    imgpaths=images;
end

for f=1:length(imgpaths)
    imgpath=imgpaths{f};
    [tmpdir,fname,ext]=fileparts(imgpath);
    if isempty(tmpdir)
        tmpdir=pwd;
    end
    outfname=[tmpdir,filesep,opt.pref,fname,ext];
    if strcmp(ext,'.gz')
      gunzip(imgpath);
      imgpath=[tmpdir,filesep,fname];
      outfname=[tmpdir,filesep,opt.pref,fname,ext];
    end
    fprintf('    --reslice image: %s --> %s\n',imgpath,outfname);
    if ~exist(outfname,'file') ||(exist(outfname,'file') && opt.overwrite==1)        
        bb=opt.bb;
        if isempty(bb) %add bondingbox optional inputs20201203
        [bb]=spm_get_bbox(imgpath);
        end
        if isnumeric(reslvox)
            if length(reslvox)==1
                vox=repmat(reslvox,[3,1]);
            else
                vox=reslvox;
            end
            matlabbatch{1}.spm.util.defs.comp{1}.idbbvox.vox = vox;
            matlabbatch{1}.spm.util.defs.comp{1}.idbbvox.bb =bb;
        elseif ischar(reslvox)
            if exist(reslvox,'file')
                [rtmpdir,rfname,rext]=fileparts(reslvox);
                if strcmp(rext,'.gz') || strcmp(rext,'.nii')
                    imgspace=reslvox;
                    [rtmpdir,rfname,rext]=fileparts(imgspace);
                    if strcmp(rext,'.gz')
                       gunzip(imgspace);
                       imgspace=[rtmpdir,filesep,rfname];
                    end
                    matlabbatch{1}.spm.util.defs.comp{1}.id.space = {imgspace};
                else
                    error('Missing image file defining the reslice space: %s\n',reslvox);
                end
            else
                error('Missing image file defining the reslice space: %s\n',reslvox);
            end
        end
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {imgpath};
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr={tmpdir};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = opt.interp;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = opt.pref;
        spm('defaults', 'FMRI');
        spm_jobman('run',matlabbatch);
        if strcmp(ext,'.gz')
          delete(imgpath);        
          gzip([tmpdir,filesep,opt.pref,fname]);
          delete([tmpdir,filesep,opt.pref,fname]);
        end
        if strcmp(ext,'.gz')
          delete(imgspace);
        end
    end
end
