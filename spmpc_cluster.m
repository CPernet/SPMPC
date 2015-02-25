function out = spmpc_cluster(what)

% computes the nb of clusters, their size, and coordinate: what == 1 (input observed data)
% computes the mean/std of the biggest clusters under H0: what == 2 (input bootstraped H0 data)

% cyril pernet 20/12/2008

global defaults;
spm_defaults

switch what

%% ------------------------------------------------------------------------
    
    case(1)
        P = spm_select(1,'image','Select the abnormal image');
        M = spm_vol(P);
        spm_check_orientations(M);
        abnormal = spm_read_vols(M);
        xmax  = M(1).dim(1);
        ymax  = M(1).dim(2);
        zmax  = M(1).dim(3);
        nbimage = size(M,1);

        t = spm_input('data threshold: % different','-1'); % threshold
        result_image = (abnormal>=t);
        nb_voxels = length(find(result_image));

        if nb_voxels>0

            % [L,NUM] = spm_bwlabel(double(result_image.*1),18); crashes
            [L,nb_clusters] = bwlabeln(result_image,18); % needs the image processing toolbox

            disp('-------------------------------------------------------')
            fprintf('%g clusters with %g%% of patients different were found',nb_clusters, t)
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

            out = [ mni_coordinates cluster_size']
            save mni_coordinates.txt out -ascii 
            
            % get clusters only
            fprintf('biggest cluster = %g voxels',max(cluster_size));
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
                Info_img = M;
                Info_img.fname = sprintf('%s%s',pwd,name);
                Info_img.description = 'thresholded classification image';
                spm_write_vol(Info_img,thresholded_image);
            end

        end % closes if nb_voxels>0


%% ------------------------------------------------------------------------

    case(2)
        P = spm_select(1,'mat','Select the H0_classification file');
        load (P)
        xmax  = size(H0_classification,2);
        ymax  = size(H0_classification,3);
        zmax  = size(H0_classification,4);
        nbimage = size(H0_classification,1);
               
        t = spm_input('data threshold: % different','-1'); % threshold
        result_image = zeros(xmax,ymax,zmax);
        index = 1;
        
        get_mask = spm_input('load mask?', +1, 'y/n');
        if get_mask == 'n'
            mask = ones(xmax,ymax,zmax);
        else
           P = spm_select(1,'image','Select mask image');
           M = spm_vol(P);
           spm_check_orientations(M);
           mask  = spm_read_vols(M);
           if M(1).dim(1) ~= xmax || M(1).dim(2) ~= ymax || M(1).dim(3) ~= zmax
               error('dimensions do not agree, error when reading the mask, error l. 105');
           end
        end
        
        
        for i=1:nbimage
            
            result_image = (squeeze(H0_classification(i,:,:,:))>=t).*mask;
            nb_voxels = length(find(result_image));
            
            if nb_voxels>0
                [L,nb_clusters] = bwlabeln(result_image,18);
                out(index) = max(length(find(L==i)));
            else
                out(index) = 0;
            end
      
            index = index+1;            
        end
    
    save max_cluster_distribution.txt out -ascii
    Alpha = 5/100;
    spm_progress_bar('Init',100,'Computing cluster inference','% boot completed')
    for l=1:1000
        [bootindex, temp] = spmpc_sample(nbimage,out);
        tempboot(l)=sum(temp)./length(temp);
        percent = (l ./ 1000) * 100;
        spm_progress_bar('Set',percent);
    end
    spm_progress_bar('Clear');
    low = round(Alpha.*l./2);
    high = l - low;
    bootsort=sort(tempboot);
    bootCI(1) = bootsort(low+1);
    bootCI(2) = bootsort(high);
    figure
    subplot(1,2,1);hist(out);title('observed max clusters');
    param1 = mle('normal',out); title2=sprintf('mean %g std %g',param1(1),param1(2));
    [a,b]=hist(out);text(b(2),max(a)-0.5,title2); grid on
    subplot(1,2,2);hist(bootsort);title('bootstrapped max clusters');
    param2 = mle('normal',bootsort); title2=sprintf('mean %g std %g',param2(1),param2(2));
    [a,b]=hist(bootsort);text(b(1),max(a)-0.5,title2); grid on
        
        
        
%% ------------------------------------------------------------------------

    otherwise


        error('specify argument 1 or 2')
end


