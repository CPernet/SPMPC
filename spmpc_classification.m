function spmpc_classification(estimate_H0, make_classification_images, P_controls, P_patients, path, alpha_value, Nboot)

% this function performs a classication of images based on CI
% - for each voxel, each scan get its' voxels at -1 0 or 1
%   percentage maps are computed and possibly created
% - the significance of the classification is tested under H0
%   H0: control and patients come from the same population (no differences)
%   Compare observed to H0 using the maximum likelihood estimate of the variance
%
% spmpc_classification(estimate_H0, make_classification_images, P_controls, path)
% estimate_H0                : 1/0 option to test or not the classication
% make_classification_images : 1/0 option to write all results as images
% P_controls                 : list control subject's images (see spm_select)
% path                       : working directory
%
% cyril pernet 02/07/2008 v2


%% initialize variables
clc
global defaults;
spm_defaults


if estimate_H0 == 1
    warndlg('the estimation of H0 is turned on, this will take time ...','long computation time expected')
    
    if alpha_value == 0 
        alpha_value = spm_input('Chosse your alpha_value level',-2,'b','1%|5%',[1,5],1);
        alpha_value = alpha_value/100;
    end
    
    if P_controls == 0
        P_controls = spm_select(Inf,'.*\.img$','Select Images the CI was computed from');
    end
    
    if Nboot == 0
        Nboot = spm_input('How Many Bootstrap were used for the CI',1);
        if Nboot==0
            disp('bye'); return

        elseif Nboot<2500
            while Nboot<2500 && Nboot~=0
                errordlg('Nboot had to be at least 2500','Nboot error')
                Nboot = spm_input('How Many Bootstraap to do (min 2500)',1);
                if Nboot==0
                    disp('bye'); return
                end
            end
        end
    end
end


%% Get the patient's data, CI, and mask

% Get patient's data
V = spm_vol(P_patients); spm_check_orientations(V);
Image = spm_read_vols(V);
xmax  = V(1).dim(1);
ymax  = V(1).dim(2);
zmax  = V(1).dim(3);
nbimage = size(V,1);

classification = zeros(5,xmax,ymax,zmax,nbimage);

% Get the confidence interval
cd(path)
try
    disp('load CI')
    load bootCI
catch
    [t,sts] = spm_select(1,'.*\.mat$','Select a 5D CI matrix');
    if sts==1
        load (t)
    else
        close all
        disp('classification aborded')
    return
    end
end

%% check dimensions

disp('check dimensions')
try
    if size(bootCI,3)~=xmax || size(bootCI,4)~=ymax || size(bootCI,5)~=zmax
    disp('error, matrix dimension must agree')
    return
    end
catch
    disp('the expected name of the CI is bootCI')
    disp('error at line 93 in spmpc_classification')
    return
end


%% Get a mask
try
    P_mask = sprintf('%s/mask.img',path)
    M = spm_vol(P_mask); disp('load mak')
    mask = spm_read_vols(M);
    mxmax  = M(1).dim(1); mymax  = M(1).dim(2); mzmax  = M(1).dim(3);
    if mxmax ~= xmax || mymax ~= ymax || mzmax ~= zmax
        disp('error, mask dimension must agree')
        return
    end
catch
    P_mask = spm_select(1,'.*\.img$','Select a mask');
    M = spm_vol(P_mask); mask = spm_read_vols(M);
    mxmax  = M(1).dim(1); mymax  = M(1).dim(2); mzmax  = M(1).dim(3);
    if mxmax ~= xmax || mymax ~= ymax || mzmax ~= zmax
        disp('error, mask dimension must agree')
        return
    end
end

%% Apply mask

disp('apply mask to images')
for n=1:nbimage
    Image(:,:,:,n) = Image(:,:,:,n) .*mask;
end


%% Perform the classification

spm_progress_bar('Init',100,'Computing','% subjects classified')
index =1;

for iteration = 1:5
    fprintf('iteration %g --------------------',iteration); disp(' ')
    down = squeeze(bootCI(iteration,1,:,:,:));
    up   = squeeze(bootCI(iteration,2,:,:,:));
    
    for s=1:nbimage
        fprintf('subject %g',s);disp(' ');
        percent=(index)/(5*nbimage)*100;
        spm_progress_bar('Set',percent); 
        
        tmp1(:,:,:,s) = squeeze(Image(:,:,:,s)) < down ;
        tmp2(:,:,:,s) = squeeze(Image(:,:,:,s)) > up ;
        tmp3 = squeeze(classification(iteration,:,:,:,s));
        tmp3(tmp1(:,:,:,s)) = -1;
        tmp3(tmp2(:,:,:,s)) = 1;
        classification(iteration,:,:,:,s) = tmp3;
        clear tmp1 tmp2 tmp3
        index = index+1;
    end 
end

disp('classication done ')


%% make intermediate percentage maps

Info_img = M;
if make_classification_images == 1

    disp('conpute percentage images ...')
    class_inferior = (classification == -1);
    class_superior = (classification == 1);
    for i=1:5
        tmp=squeeze(class_inferior(i,:,:,:,:));
        inf(i,:,:,:) = squeeze((sum(tmp,4).*100)./nbimage);
        tmp=squeeze(class_superior(i,:,:,:,:));
        sup(i,:,:,:) = squeeze((sum(tmp,4).*100)./nbimage);
    end
    inferior = squeeze(sum(inf,1)./5).*mask;
    superior = squeeze(sum(sup,1)./5).*mask;
    clear tmp inf sup class_inferior class_superior

    normal = 100-inferior-superior;
    tmp = (normal == 100); % to get 0 outside the brain
    normal(tmp) = 0;
    normal = normal.*mask;

    name = '/normal.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.description = '% subject within CI';
    spm_write_vol(Info_img,normal);
    name = '/superior.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.description = '% subject sup to the CI';
    spm_write_vol(Info_img,superior);
    name = '/inferior.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.description = '% subject inf to the CI';
    spm_write_vol(Info_img,inferior);
end

clear normal inferior superior tmp
spm_progress_bar('Clear')


%% final classication step
% need to cut the 5D matrix for memory reasons ..
% needed have to put that in a loop according to the 
% number of CI we want to have + some computers will crash
% to heavy need to loop in z at least

disp('computing the final average classification')

new_classification1 = squeeze(classification(1,:,:,:,:));
new_classification1(squeeze(classification(1,:,:,:,:) == -1))=1;
abnormal1 = (sum(new_classification1,4)./nbimage*100);
clear new_classification1
disp('classifcation 1 done')

new_classification2 = squeeze(classification(2,:,:,:,:));
new_classification2(squeeze(classification(2,:,:,:,:) == -1))=1;
abnormal2 = (sum(new_classification2,4)./nbimage*100);
clear new_classification2
disp('classifcation 2 done')

new_classification3 = squeeze(classification(3,:,:,:,:));
new_classification3(squeeze(classification(3,:,:,:,:) == -1))=1;
abnormal3 = (sum(new_classification3,4)./nbimage*100);
clear new_classification3
disp('classifcation 3 done')

new_classification4 = squeeze(classification(4,:,:,:,:));
new_classification4(squeeze(classification(4,:,:,:,:) == -1))=1;
abnormal4 = (sum(new_classification4,4)./nbimage*100);
clear new_classification4
disp('classifcation 4 done')

new_classification5 = squeeze(classification(5,:,:,:,:));
new_classification5(squeeze(classification(5,:,:,:,:) == -1))=1;
abnormal5 = (sum(new_classification5,4)./nbimage*100);
clear new_classification5
disp('classifcation 5 done')

abnormal = (abnormal1+abnormal2+abnormal3+abnormal4+abnormal5)./5;
name = '/abnormal.img';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.descrip = '% abnormal';
spm_write_vol(Info_img,abnormal);

if make_classification_images == 1
    name = '/abnormal_1.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.descrip = '% abnormal 1st CI';
    spm_write_vol(Info_img,abnormal1);
        
    name = '/abnormal_2.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.descrip = '% abnormal 1st CI';
    spm_write_vol(Info_img,abnormal2);
    
    name = '/abnormal_3.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.descrip = '% abnormal 1st CI';
    spm_write_vol(Info_img,abnormal3);
    
    name = '/abnormal_4.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.descrip = '% abnormal 1st CI';
    spm_write_vol(Info_img,abnormal4);
    
    name = '/abnormal_5.img';
    Info_img.fname = sprintf('%s%s',path,name);
    Info_img.descrip = '% abnormal 1st CI';
    spm_write_vol(Info_img,abnormal5);
end

clear abnormal1 abnormal2 abnormal3 abnormal4 abnormal5
classification = squeeze(sum(classification,1)./5);
save classification classification
disp('classifcation finished ')


%% find areas with X% abnormal

test_threshold = 'y';

while test_threshold == 'y'
    t = spm_input('data threshold: % different','-1'); % threshold
    result_image = (abnormal>=t);
    nb_voxels = length(find(result_image));

    if nb_voxels>0

        % [L,NUM] = spm_bwlabel(double(result_image.*1),18); crashes
        [L,nb_clusters] = bwlabeln(result_image,18); % needs the image processing toolbox

        disp('-------------------------------------------------------')
        if t == 100
            fprintf('%g clusters with 100%% of patients different were found',nb_clusters)
        else
            fprintf('%g clusters with at least %g%% of patients different were found',nb_clusters, t)
        end
        disp(' ')
        disp('-------------------------------------------------------')

        % voxel_coordinates = zeros(nb_voxels,3);
        mni_coordinates = zeros(nb_clusters,3);
        index = 0;
        for i=1:nb_clusters
            tmp = (L==i);
            for z=1:zmax
                [c1,c2]=find(tmp(:,:,z));
                if isempty(c1) == 0
                    for n=1:length(c1)
                        index = index+1;
                        voxel_coordinates(index,:)= [c1(n),c2(n),z];
                    end
                end
            end
            index = 0;
            cluster_size(i) = length(find(L==i));
            mni_coordinates(i,:)= (M.mat(1:3,:)* [mean(voxel_coordinates,1) 1]')';
               
        end

        coord = [ mni_coordinates cluster_size']
        save mni_coordinates.txt coord -ascii
        fprintf('maximum cluster size %g',max(cluster_size)); disp(' ')

        % get clusters only
        choice = spm_input('write tresholded/cluster image?','-1','y/n');

        if choice =='y'

            cluster_threshold = spm_input('cluster size?',1);
            for i=1:nb_clusters
                if cluster_size(i)<(cluster_threshold-1)
                    L(L==i)=0;
                end
            end
            L(L~=0)=1;
            thresholded_image = abnormal.*L;
            name = '/thresholded_image.img';
            Info_img.fname = sprintf('%s%s',pwd,name);
            Info_img.description = 'thresholded classification image';
            spm_write_vol(Info_img,thresholded_image);
        end


        % compute correlations across voxels
        choice = spm_input('compute correlations across voxels?','-1','y/n');

        if choice =='y'
            index = 1; spm_progress_bar('Init',100,'Getting original values','% voxels')

            for z=1:zmax-2
                percent = z/(zmax-2)*100;
                [x,y]=find(abnormal(:,:,z)>=t);
                test = sum(sum(isempty(x)));

                if test == 0
                    for q=1:length(x)
                        location{index} = [x(q),y(q),z]; % x,y,z coordinate
                        distribution(:,index) = squeeze(Image(x(q),y(q),z,:)); % original values
                        remaining_vectors(:,index) = squeeze(classification(x(q),y(q),z,:)); % classification
                        index=index+1;
                    end
                end
                spm_progress_bar('Set',percent);
            end
            spm_progress_bar('Clear')

            % have a look at how voxels correlates
            figure
            subplot(1,2,1);imagesc(corr(remaining_vectors,'type','spearman')); title('Spearman correlations (class)')
            subplot(1,2,2);imagesc(corr(distribution)); title('Pearson correlations (distrib)')
        end
        
    else
        
        fprintf('0 voxel superior or equal to the %g%%');

    end % closes if nb_voxels>0

    test_threshold = spm_input('recompute for another threshold?',-1,'y/n');
end

%% testing H0

if estimate_H0 == 1
    disp('start testing the null hypothesis')
    disp('compute 100 H0 samples')

    save testH0 P_controls P_patients P_mask xmax ymax zmax alpha_value Nboot
    clear all; load testH0
    parameters = spmpc_testH0(P_controls, P_patients, P_mask, xmax, ymax, zmax, alpha_value, Nboot);

    disp('compute the p values ..')

    name = sprintf('%s/abnormal.img',pwd)
    V = spm_vol(name); spm_check_orientations(V);
    observed = spm_read_vols(V);
    probability = zeros(xmax,ymax,zmax);
    spm_progress_bar('Init',100,'estimating normal parameters','% slices completed')

    for k=1:zmax
        percent = k/zmax*100

        for j= 1:ymax;
            for i= 1:xmax;

                if observed(i,j,k) > 0
                    probability(i,j,k) = pdf('Normal',observed(i,j,k),parameters(i,j,k,1),parameters(i,j,k,2));
                end

            end % close the loop for x
        end % close the loop for yj

        spm_progress_bar('Set',percent);

    end % close the loop for z

    
    Info_img = V;
    Info_img.fname = sprintf('%s/probability_map.img',pwd);
    Info_img.descrip = 'prob of being abnormal';
    spm_write_vol(Info_img,probability);

end

    