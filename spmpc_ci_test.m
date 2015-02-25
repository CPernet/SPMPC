function spmpc_ci_test(path,P)

% this function evaluates as many data points as you want
% by default 5000 bootstraap are performed
% this gives a flavor of the stability of the CI
%
% cyril pernet 27/06/08 v2


%% initialize variables
global defaults
spm_defaults
Nboot = 5000;

%% get data
if P == 0
    [P,sts] = spm_select(Inf,'.*\.img$','Select Images to compute CI from');
end

while isempty(P) == 1 && sts == 1
    [P,sts] = spm_select(Inf,'.*\.img$','Select Images to compute CI from');
end

if sts == 1 

    spm_orientations(P); V = spm_vol(P); spm_check_orientations(V);
    Image = spm_read_vols(V); nbimage = size(V,1);

    %% compute mean, std, etc ..
    choice = spm_input('load summary statistic images from the folder?','-1','y/n');
    if choice =='n'
        binary_img = Image;
        disp('computing the mask');
        for n=1:nbimage % for each image
            for z = 1:V(1).dim(3) % for each slice
                [N,X]=hist(binary_img(:,:,z));
                binary_img(:,:,z,n)=Image(:,:,z,n) >= X(2);
            end
        end

        mask = sum(binary_img,4);
        mask = (mask == nbimage);

        disp('computing the mean image');
        mean_image = zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3)); mean_image = mean(Image,4) .* mask ;
        disp('computing the std image');
        std_image = zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3)); std_image = std(Image,0,4) .* mask ;
        disp('computing the kurtosis image');
        kurtosis_image = zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3)); kurtosis_image = kurtosis(Image,0,4) .* mask;
        disp('computing the skewness image');
        skewness_image = zeros(V(1).dim(1),V(1).dim(2),V(1).dim(3)); skewness_image = skewness(Image,0,4) .* mask;

        disp('write images');
        Info_img = V(1);
        name = '/mean_image.img';
        Info_img.fname = sprintf('%s%s',path,name);
        Info_img.descrip = 'mean_image';
        spm_write_vol(Info_img,mean_image);

        name = '/std_image.img';
        Info_img.fname = sprintf('%s%s',path,name);
        Info_img.descrip = 'std_image';
        spm_write_vol(Info_img,std_image);

        name = '/kurtosis_image.img';
        Info_img.fname = sprintf('%s%s',path,name);
        Info_img.descrip = 'kurtosis_image';
        spm_write_vol(Info_img,kurtosis_image);

        name = '/skewness_image.img';
        Info_img.fname = sprintf('%s%s',path,name);
        Info_img.descrip = 'skewness image';
        spm_write_vol(Info_img,skewness_image);
    end

    %% do the bootstrap for a given region
    stop = 0;
    while stop == 0
        choice = spm_input('mean, std, skew, kurt? 1->4','-1');
        if choice ==1
            name = '/mean_image.img';
            Info_img.fname = sprintf('%s%s',path,name);
        elseif choice ==2
            name = '/std_image.img';
            Info_img.fname = sprintf('%s%s',path,name);
        elseif choice ==3
            name = '/skewness_image.img';
            Info_img.fname = sprintf('%s%s',path,name);
        elseif choice ==4
            name = '/kurtosis_image.img';
            Info_img.fname = sprintf('%s%s',path,name);
        else
            disp('bye')
            stop = 1;
            close all
            return
        end

        close
        spm_image('init',Info_img.fname)

        coord = spm_input('enter a voxel coordinates',1);
        data = Image(coord(1),coord(2),coord(3),:);

        % compute the bootstrap on Y
        Alpha = spm_input('Chosse your alpha_value level',-2,'b','1%|5%',[1,5],1);
        Alpha = Alpha/100;
        index = 1;
        spm_progress_bar('Init',100,'Computing CI','% boot completed')

        figure('Name','evolution of the bootstrap')
        bootCI = zeros(2,Nboot);
        tempboot(1:Nboot) = 0;

        for l=1:Nboot
            fprintf('boot %g',l);
            disp(' ');
            [bootindex, temp] = spmpc_sample(nbimage,data);
            tempboot(l)=sum(temp)./length(temp);
            percent = (l ./ Nboot) * 100;
            spm_progress_bar('Set',percent);
            low = round(Alpha.*l./2);
            high = l - low;
            bootsort=sort(tempboot(1:index));
            bootCI(1,index) = bootsort(low+1);
            bootCI(2,index) = bootsort(high);
            if l<Nboot
                subplot(1,2,1);
                plot( bootCI(2,1:index) - bootCI(1,1:index) )
                subplot(1,2,2)
                hist(tempboot(1:index))
            elseif l==Nboot
                subplot(2,2,1)
                hist(data, length(data));
                title('data'); grid on

                subplot(2,2,2)
                plot(bootCI(1,:));
                hold on; plot(bootCI(2,:));
                [h,p,ci]=ttest(squeeze(data));
                plot(repmat(ci(1),5000,1),'r')
                plot(repmat(ci(2),5000,1),'r')
                title('evolution of bootstrap CIs');
                grid on

                subplot(2,2,3);
                plot( bootCI(2,1:index) - bootCI(1,1:index) )
                title('CI size'); grid on
                subplot(2,2,4)
                hist(tempboot(1:index))
                title('histogram of the means');
                grid on
            end

            drawnow
            index=index+1;
        end

        spm_progress_bar('Clear');

        p = spm_input('do you want to test another ponit?','-1','y/n');
        if p=='n'
            disp('bye')
            close Figure 2
            stop = 1;
            return
        else
            try
                close Figure 1
                close Figure 2
            catch
                close Figure 2
            end
            clear tempboot data
        end
    end

end
