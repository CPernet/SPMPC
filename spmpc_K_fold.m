function spmpc_K_fold(path,P_controls,P_patients,alpha_value,Nboot)

% this function is to compute a 3 fold cross validation
% control data are split in 3, 2/3 are used in turn to compute the CI
% then patients are classified - after the 3 rounds, results are averaged
%
% cyril pernet 05/11/2008 v4

%% --------------
%   get the data


P_controls=0;
P_patients=0;
alpha_value=0;
Nboot=0;
one_sample_ci = 0;
ci_images = 0;


global defaults
spm_defaults;


% control subjects
if P_controls == 0
    P_controls= spm_select(Inf,'.*\.img$','Select control images');
end
V1 = spm_vol(P_controls);
spm_check_orientations(V1);
xmax  = V1(1).dim(1);
ymax  = V1(1).dim(2);
zmax  = V1(1).dim(3);
nbimage_controls = size(V1,1);

% patients
if P_patients == 0
    P_patients= spm_select(Inf,'.*\.img$','Select dyslexics images');
end
V2 = spm_vol(P_patients);
spm_check_orientations(V2);
nbimage_dyslexics = size(V2,1);

% p level
if alpha_value == 0
    alpha_value = spm_input('Chosse your alpha_value level',-1,'b','1%|5%',[1,5],1);
    alpha_value = alpha_value/100;
end

% nb of bootstraps
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


%% create the K-kold directory

mkdir ('3_fold_cross_validation')
cd ('3_fold_cross_validation')

%% ---------------
% get a mask

p = spm_input('do you want to load a mask ?','-1','y/n');

if p =='n' % compute a mask

    nbimage = nbimage_controls + nbimage_dyslexics;
    Controls_Images = spm_read_vols(V1);
    Dyslexics_Images = spm_read_vols(V2);
    Image = zeros(size(Controls_Images,1),size(Controls_Images,2),size(Controls_Images,3),nbimage);
    Image(:,:,:,1:size(Controls_Images,4)) = Controls_Images;
    Image(:,:,:,size(Controls_Images,4)+1:nbimage) = Dyslexics_Images;

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
    Info_img = V1(1);
    name = '/mask.img';
    Info_img.fname = sprintf('%s%s',pwd,name);
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



%% 3 fold - create random groups

K = floor(nbimage_controls / 3);
a = ceil( rand(nbimage_controls,1)*nbimage_controls );
b = ceil( rand(nbimage_controls,1)*nbimage_controls );
selection = [a b];
selection = selection(selection ~= 0);
index = 1;
gp12 = zeros(2*K,1);
stop =0;
while stop ==0
    for i=1:length(selection)
        if selection(i)<=nbimage_controls && any(gp12==selection(i)) == 0 && selection(i) ~= 0
            gp12(index) = selection(i);
            index = index +1;
            if index == 2*K
                stop =1;
            end
        end
    end
end

gp1 = gp12(1:K);
gp2 = gp12(K+1:2*K);
gp12 = [gp1 ;gp2];

gp3 = zeros(K,1);
index = 1;
for i=1:nbimage_controls
    test = sum(i == gp12);
    if test == 0
        gp3(index) = i;
        index = index + 1;
    end
end


%% conpute the K folds


for i= 1:3

    % create a directory
    name = sprintf('fold_%g',i);
    mkdir(name)
    cd(name)

    if i== 1

        % make the control gp and get the CI
        P = [P_controls(gp1,:); P_controls(gp2,:)];
        [Alpha, Nboot] = spmpc_ci(one_sample_ci,ci_images,P,pwd,alpha_value,Nboot,mask);
        % take the remaining controls and classify
        P = [P_controls(gp3,:)];
        path = pwd; classify_controls(P, path);
        % classify patients
        path = pwd; classify_patients(P_patients, path);
        cd ..

    elseif i == 2
        P = [P_controls(gp1,:); P_controls(gp3,:)];
        [Alpha, Nboot] = spmpc_ci(one_sample_ci,ci_images,P,pwd,alpha_value,Nboot,mask);
        % take the remaining controls and classify
        P = [P_controls(gp2,:)];
        path = pwd; classify_controls(P, path);
        % classify patients
        path = pwd; classify_patients(P_patients, path);
        cd ..

    elseif i ==3
        P = [P_controls(gp2,:); P_controls(gp3,:)];
        [Alpha, Nboot] = spmpc_ci(one_sample_ci,ci_images,P,pwd,alpha_value,Nboot,mask);
        % take the remaining controls and classify
        P = [P_controls(gp1,:)];
        path = pwd; classify_controls(P, path);
        % classify patients
        path = pwd; classify_patients(P_patients, path);
        cd ..

    end

end


%% ----------------------
function classify_controls(P, path);

V = spm_vol(P); spm_check_orientations(V);
Info_img = V; Image = spm_read_vols(V);
xmax  = V(1).dim(1); ymax  = V(1).dim(2);
zmax  = V(1).dim(3); nbimage = size(V,1);
disp('apply mask to images')
for n=1:nbimage
    Image(:,:,:,n) = Image(:,:,:,n) .*mask;
end
classification = zeros(5,xmax,ymax,zmax,nbimage);
spm_progress_bar('Init',100,'Computing','% subjects classified FOLD 1')
index =1; load bootCI

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

clear bootCI
disp('computing the average classification')

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

mkdir('controls'); cd('controls'); path = pwd;
abnormal = (abnormal1+abnormal2+abnormal3+abnormal4+abnormal5)./5;
name = '/abnormal_controls.img';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.descrip = '% abnormal_controls';

spm_write_vol(Info_img,abnormal);
clear abnormal1 abnormal2 abnormal3 abnormal4 abnormal5
classification = squeeze(sum(classification,1)./5);
save classification_controls classification
cd ..

%% ---------------------------------------------
function classify_patients(P_patients, path);


V = spm_vol(P_patients); spm_check_orientations(V);
Info_img = V; Image = spm_read_vols(V);
xmax  = V(1).dim(1); ymax  = V(1).dim(2);
zmax  = V(1).dim(3); nbimage = size(V,1);
disp('apply mask to images')
for n=1:nbimage
    Image(:,:,:,n) = Image(:,:,:,n) .*mask;
end
classification = zeros(5,xmax,ymax,zmax,nbimage);
spm_progress_bar('Init',100,'Computing','% subjects classified FOLD 1')
index =1; load bootCI

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

clear bootCI
disp('computing the average classification')

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

mkdir('patients');cd('patients'); path = pwd;
abnormal = (abnormal1+abnormal2+abnormal3+abnormal4+abnormal5)./5;
name = '/abnormal_patients.img';
Info_img.fname = sprintf('%s%s',path,name);
Info_img.descrip = '% abnormal';
spm_write_vol(Info_img,abnormal);
clear abnormal1 abnormal2 abnormal3 abnormal4 abnormal5
classification = squeeze(sum(classification,1)./5);
save classification classification
cd ..
