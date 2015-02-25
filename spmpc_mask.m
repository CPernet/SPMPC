function spmpc_mask

% silly function to create a mask
% cyril pernet 26/06/08

global defaults
spm_defaults
P= spm_select(Inf,'.*\.img$','Select Images to compute the mask from');
spm_orientations(P);
V = spm_vol(P);
spm_check_orientations(V);
Image = spm_read_vols(V);
xmax  = V(1).dim(1);
ymax  = V(1).dim(2);
zmax  = V(1).dim(3);
nbimage = size(V,1);


pp = spm_input('threshold implicit masking ?','-1','y/n');

if pp == 'n' % threshold
    threshold = spm_input('threshold value ?','-1');

    for n=1:nbimage
        binary_img(:,:,:,n) = Image(:,:,:,n) > threshold;
    end

else % assume threshold as the 2st bin of the histogramm

    for n=1:nbimage
        binary_img(:,:,:,n)=Image(:,:,:,n);
        for z = 1:zmax
            [N,X]=hist(binary_img(:,:,z));
            binary_img(:,:,z,n)=Image(:,:,z,n) > X(2);
        end
    end

end


mask = sum(binary_img,4);
mask = (mask == nbimage);

% save the mask
Info_img = V(1);
path = uigetdir(pwd,'select saving directory');
name = '/mask.img';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.descrip = 'mask image';
spm_write_vol(Info_img,mask);
