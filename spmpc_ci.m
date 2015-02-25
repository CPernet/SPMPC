function [Alpha, Nboot] = spmpc_ci(one_sample_ci,ci_images,P,path,Alpha,Nboot,mask)

% this function create a boostraped confidence interval from a single group
% please note that the more bootstrap you do the more acurate is the CI but
% also the more demanding it is for your computer; this can take several
% hours - the function returns a 5D matrix boorCI(5,2,x,y,z)
% 5 CI for the last 2000 resamples
% 2 for lowe and upper bounds
% x,y,z spatial coordinates
%
% one_sample_ci : 0/1 binary decision if one wants to compute the one sample t test CI
% ci_images     : 0/1 binary decision if one wants to generate image(s) of the CI(s) size(s)
% P             : list the images to process (see P = spm_select)
% path          : indicate the working directory - the bootCI.mat 5D matrix will be saved there
% alpha_value         : returns the alpha_value level used to compute the CI, this is used when evaluating H0 (spmpc_testH0)
% 
% Cyril Pernet 30/06/2008 v3


%% inintialize variables
global defaults
spm_defaults

%%   get the data
% ---------------

cd (path)
spm_orientations(P);
V = spm_vol(P);
spm_check_orientations(V);
Image = spm_read_vols(V);
xmax  = V(1).dim(1);
ymax  = V(1).dim(2);
zmax  = V(1).dim(3);
nbimage = size(V,1);

bootCI=zeros(5,2,xmax,ymax,zmax);
if one_sample_ci == 1; CI=zeros(2,xmax,ymax,zmax); end


%% get a mask to reduce the nb of computations
% --------------------------------------------

if mask == 0

    p = spm_input('do you want to load a mask ?','-1','y/n');

    if p =='n' % compute a mask

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
                    [N,X]=hist(binary_img(:,:,z)); % I wish I could vectorize this ..
                    binary_img(:,:,z,n)=Image(:,:,z,n) > X(2);
                end
            end

        end


        mask = sum(binary_img,4);
        mask = (mask == nbimage);

        % save the mask
        Info_img = V(1);
        name = '/mask.img';
        Info_img.fname = sprintf('%s%s',newpath,name);
        Info_img.descrip = 'mask image';
        spm_write_vol(Info_img,mask);

    else % load a mask

        [M, sts] = spm_select(1,'.*\.img$','Select a mask');
        if sts == 1
            M = spm_vol(M);
            mask = spm_read_vols(M);
            mxmax  = M(1).dim(1);
            mymax  = M(1).dim(2);
            mzmax  = M(1).dim(3);
            if mxmax ~= xmax || mymax ~= ymax || mzmax ~= zmax
                disp('error, mask dimension must agree')
                return
            end
        else
            return
        end

    end
end


%% bootstrap parameters
% ---------------------

if Nboot == 0
    Nboot = spm_input('How Many Bootstrap (min 2500)',1);
    if Nboot==0
        disp('bye'); return

    elseif Nboot<2500
        while Nboot<2500 && Nboot~=0
            errordlg('Nboot has to be at least 2500','Nboot error')
            Nboot = spm_input('How Many Bootstraap to do (min 2500)',1);
            if Nboot==0
                disp('bye'); return
            end
        end
    end
end

Nboot = abs(Nboot); % in case of typo error

if Alpha == 0
    Alpha = spm_input('Chosse your alpha_value level',-2,'b','1%|5%',[1,5],1);
    Alpha = Alpha/100;
end
threasholds = [Nboot-2000; Nboot-1500; Nboot-1000; Nboot-500; Nboot];
low = round(Alpha.*threasholds./2);
high = threasholds - low;

% get a bootstraap index used for all voxels at each bootstrap

for l=1:Nboot
    bootindex(:,l) = ceil( rand(nbimage,1)*nbimage ); % same problem as bellow if l=1:Nboot, bootindex(:,l) isn't working
end


%% get bootstap CI
% ----------------

spm_progress_bar('Init',100,'Computing bootstrapped CI','% slices completed')
l=1:Nboot;

for k= 1:zmax;

    fprintf('slice %g',k); disp(' '); percent = k/zmax*100;

    for j= 1:ymax;
        for i= 1:xmax;

            if mask(i,j,k) > 0
                tempdata = squeeze(Image(i,j,k,bootindex(:,l))); % vector of bootstrapped data - I can't it as matrix
                data = zeros(nbimage,Nboot); % matrix for bootstrapped data to fill
                data(:,1) = tempdata(1:nbimage); % reorganiz data - maybe there is a clever way ?
                index = 1;
                for c = 2:Nboot
                    data(:,c) = tempdata( (nbimage*(c-1))+1 : nbimage*c);
                end
                tempboot = mean(data);
                bootsort=sort(tempboot);
                bootCI(:,1,i,j,k) = bootsort(low+1);
                bootCI(:,2,i,j,k) = bootsort(high);
                clear tempboot
            end

        end % close the loop for x
    end % close the loop for y

    spm_progress_bar('Set',percent);

end % close the loop for z

spm_progress_bar('Clear');
disp('bootstrapped CI computation done')
save bootCI bootCI


%% do the same analysis for the t-test ci
% ---------------------------------------

if one_sample_ci == 1
    warning off
    j=1:ymax; k=1:xmax;
    h(k,j)=0; p(k,j)=0;
    spm_progress_bar('Init',100,'Computing t-test CI','% slices completed')
    disp(' compute basic CI');
    for i=1:zmax
        fprintf('slice %g on %g',i,zmax); disp(' ');percent = i/zmax*100;
        data = squeeze(Image(k,j,i,:));
        for n=1:nbimage;newdata(n,:,:)=data(:,:,n);end % reshape - is there a better way?
        [h(k,j),p(k,j),CI(:,k,j,i)] = ttest(newdata,0,Alpha,'both');
        CI(1,k,j,i) = squeeze(CI(1,k,j,i)).*mask(k,j,i); % clear the CI for 0s data - back to 0s
        CI(2,k,j,i) = squeeze(CI(2,k,j,i)).*mask(k,j,i);
        spm_progress_bar('Set',percent);
    end % close the loop for z

    spm_progress_bar('Clear');
    disp('one sample t-test CI computation done')
    
    save ttestCI CI
    warning on

end


%% generate images
% ----------------

if ci_images == 1
    lowboot  = squeeze(min(squeeze(bootCI(:,1,:,:,:)))); % the min in x,y,z
    highboot = squeeze(max(squeeze(bootCI(:,2,:,:,:)))); % the max in x,y,z
    bootCI_size = highboot - lowboot;
    Info_img = V(1);
    name = '/bootCI_size.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.descrip = 'size of the boot CI';
    spm_write_vol(Info_img,bootCI_size);

    if one_sample_ci == 1
        CI_size(:,:,:) = squeeze(CI(2,:,:,:)) - squeeze(CI(1,:,:,:));
        diff(:,:,:) = CI_size - bootCI_size;
        name = '/ttest_CI_size.img';
        Info_img.fname = sprintf('%s%s',path,name);
        Info_img.descrip = 'size of the one sample ttest CI';
        spm_write_vol(Info_img,CI_size);
        name = '/CI_comparison.img';
        Info_img.fname = sprintf('%s%s',path,name);
        Info_img.descrip = 'difference between bootCI and ttest CI';
        spm_write_vol(Info_img,diff);
    end
end
