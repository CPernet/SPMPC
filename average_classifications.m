function average_classifications

clear all % need some memory here
global defaults
spm_defaults
load image_parameters
H0_classification = zeros(Nboot_sample,xmax,ymax,zmax);

disp(' compute the final classications ....')
disp(' averaging across CI ................')
disp(' ')

for data=1:Nboot_sample
    matrix = sprintf('classification%g.mat',data);
    fprintf('processing %s',matrix); disp(' ');
    load(matrix)
    
    new_classification1 = squeeze(classification(1,:,:,:,:)); % first CI = 4D matrix
    new_classification1(squeeze(classification(1,:,:,:,:) == -1))=1; % make it binary 0/1
    abnormal1 = (sum(new_classification1,4)./size(classification,5)*100); % sum across subjects
    clear new_classification1

    new_classification2 = squeeze(classification(2,:,:,:,:));
    new_classification2(squeeze(classification(2,:,:,:,:) == -1))=1;
    abnormal2 = (sum(new_classification2,4)./size(classification,5)*100);
    clear new_classification2

    new_classification3 = squeeze(classification(3,:,:,:,:));
    new_classification3(squeeze(classification(3,:,:,:,:) == -1))=1;
    abnormal3 = (sum(new_classification3,4)./size(classification,5)*100);
    clear new_classification3

    new_classification4 = squeeze(classification(4,:,:,:,:));
    new_classification4(squeeze(classification(4,:,:,:,:) == -1))=1;
    abnormal4 = (sum(new_classification4,4)./size(classification,5)*100);
    clear new_classification4

    new_classification5 = squeeze(classification(5,:,:,:,:));
    new_classification5(squeeze(classification(5,:,:,:,:) == -1))=1;
    abnormal5 = (sum(new_classification5,4)./size(classification,5)*100);
    clear new_classification5

    H0_classification(data,:,:,:) = (abnormal1+abnormal2+abnormal3+abnormal4+abnormal5)./5;
    clear abnormal1 abnormal2 abnormal3 abnormal4 abnormal5 classification matrix
    % delete(matrix)
end
    
save new_H0_classification H0_classification

