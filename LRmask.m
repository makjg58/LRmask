function LRmask(fnms)

if ~exist('fnms','var')
 	%fnms = spm_select(inf,'image','Select image[s] for NaN removal'); 
    fnms = spm_select(inf,'^.*\.(gz|voi|img|nii)$','Select images to LR mask');
end
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    [fnm, isGz] = unGzSub (fnm);
    findMidline(fnm);
    midlineMask(fnm);
    if isGz %delete uncompressed image: FSL does not allow 'img.nii' and 'img.nii.gz' to exist in same folder
        delete(fnm);
    end
end


function findMidline(fnm)
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
[pth, nam, ext] = spm_fileparts(hdr.fname);
fname_flip = fullfile(pth, ['LR', nam, ext]);
hdr_flip = hdr;
hdr_flip.fname = fname_flip;
hdr_flip.mat = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * hdr_flip.mat;
spm_write_vol(hdr_flip,img); 
%coregister data
hdr_flip = spm_vol(fname_flip); 
x  = spm_coreg(hdr, hdr_flip); 
%apply half of transform to find midline
x  = (x/2); 
M = spm_matrix(x);
MM = spm_get_space(hdr.fname);
spm_get_space(hdr.fname, M*MM); %reorient flip
delete(fname_flip);
%end doFlip()

function midlineMask(fnm)
[pth, nm, ext] = spm_fileparts(fnm);
prefix = ['a', 'b'];
rng('default')
rng('shuffle')
r = randperm(2)
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
%nul     = [0 -1.1 0.98 120]; %right slices
nul = [1 0 0 0];
d       = [size(img) 1];
[i,j,k] = ndgrid(1:d(1),1:d(2),1:d(3));
nul1    = nul(1,:)*hdr.mat;
%save Left masked
msk     = nul1(1)*i + nul1(2)*j + nul1(3)*k + nul1(4) < 0;
imgL = img .* msk;
hdr.fname = fullfile(pth, [prefix(r(1)) nm ext]); 
hdr.descrip = 'L';
spm_write_vol(hdr,imgL);
%save Right masked
msk     = nul1(1)*i + nul1(2)*j + nul1(3)*k + nul1(4) > 0;
imgR = img .* msk;
imgR = flip(imgR,1);
hdr.descrip = 'Rflipped';
hdr.fname = fullfile(pth, [prefix(r(2)) nm ext]);  
spm_write_vol(hdr,imgR);
%midlineMask

function [fnm, isGz] = unGzSub (fnm)
isGz = false;
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz -> .nii
    isGz = true;
    fnm = char(gunzip(fnm));    
elseif strcmpi(ext,'.voi') %.voi -> .nii
    isGz = true;
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end;  
%end unGzSub()
