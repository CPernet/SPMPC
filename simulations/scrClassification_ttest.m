% In this script the comparison between a classification and a t-test with
% normally distributed univariate samples will be compared. Concretely,
% Control and Test sets sampled from a Normal distribution with given  mean 
% and std (in the Parameters section) will be used. Concretely:
%
% a) Classification: The control set will be used to create a bootstrap CI
% of the means. After that, the classification will consist on looking if
% the mean of each test set falls into the CI or not.
%
% b) t-test: Will assess whether the mean of the test sets can be assumed 
% to be equal to mu or not.
%

clear all
clc

%% -- Parameters --

% mu_control = 130;
% mu_pat = 150;
mu = 130;
sigma = 20;
delta = 0:0.5:10;
numDeltas = length(delta);

numPatientSets = 200;
numSampControl = 1000;
numSampPatients = 1000;
nBoot = 5000;
alpha = 0.05;

resultsClassif = zeros(1, numDeltas);
resultsTtest = zeros(1, numDeltas);

%% -- Generate control and test samples -- 

% The control and tests are samples from a Normal distribution with mean mu
% (and mu+alpha in the case of the tests) and std sigma. First, the control 
% group is generated, and the bootstrap CI is computed. Then, the whole 
% process of generating the control sets, and classification / t-test is 
% done once per delta, storint the results in the vectors resultsClassif 
% and resultsTtest

%Sample control set
xControl = normrnd(mu,sigma, [1, numSampControl]);

% 1 get the sorted bootstrapped means
bootSort = sort(mean(xControl(randi(numSampControl,numSampControl,nBoot)),1));
% 2 CI
low = round(alpha*nBoot/2);
high = nBoot-low;
ci = [bootSort(low) bootSort(high)];

% another option is to use maximum likelihood to estimate the underlying distribution 
% and use likelihood intervals ; we can bootstrap to have better estimates
% ML is the best solution to estime a parameter
phat = mle(xControl);
% most data are at */- 2sigma that is
ci_ml = [phat(1)-2*phat(2) phat(1)+2*phat(2)];

% a boostrapped version could be ??
phatb = NaN(nBoot,2);
for b=1:nBoot
phatb(b,:) = mle(xControl(randi(numSampControl,numSampControl,1)));
end
ci_mlb = sort([phatb(:,1)-2*phatb(:,2) phatb(:,1)+2*phatb(:,2)]);
ci_mlb = [ci_mlb(low,1) ci_mlb(high,2)];

    nb_abnormal_controls_ci = (sum(xControl<ci(1))+sum(xControl>ci(2))) / numSampControl * 100;
    nb_abnormal_controls_ml = (sum(xControl<ci_mlb(1))+sum(xControl>ci_mlb(2))) / numSampControl * 100;

for i=1:numDeltas
    fprintf('running Delta = %g\n',delta(i));
    xTest =  normrnd(mu+delta(i), sigma, [numPatientSets, numSampPatients]);
    % testMeans = mean(xTest, 2);
    
    %% -- Classification by bootstrap CI --
    
    %     yClassif = zeros(size(testMeans));
    %     yClassif(testMeans<ci(1)) = -1;
    %     yClassif((ci(1)<=testMeans) & (testMeans<=ci(2))) = 0;
    %     yClassif(ci(2)<testMeans) = 1;
    %
    %     % Compute proportion of rejects
    %     resultsClassif(i) = 100*sum(yClassif~=0)/length(yClassif);
    
    % Cyril - the classification run per subject, we are not comparing means
    nb_abnormal_patients_ci{i} = sum(xTest<ci(1),1)+sum(xTest>ci(2),1);
    nb_abnormal_patients_ml{i} = sum(xTest<ci_mlb(1),1)+sum(xTest>ci_mlb(2),1);
    
    % now we use bootstrap to evaluate if this number of abnormal patients
    % is signifant (the null hypothesis is controls = patients so we simply
    % do many time the same operation but centering the data)
    H0_xControl = xControl-mean(xControl);
    H0_xTest = xTest - repmat(mean(xTest,1),200,1);
    for b= 1:600
        % compute a new ci by resampling H0 controls
        tmp = H0_xControl(randi(numSampControl,numSampControl,nBoot));
        
        bootSort = sort(mean(tmp,1));
        ci = [bootSort(low) bootSort(high)];
        
        phatb = NaN(nBoot,2); 
        parfor bb=1:nBoot
            phatb(bb,:) = mle(tmp(:,bb));
        end
        ci_mlb = sort([phatb(:,1)-2*phatb(:,2) phatb(:,1)+2*phatb(:,2)]);
        ci_mlb = [ci_mlb(low,1) ci_mlb(high,2)];
        
        % compute a new classification by resampling H0 patients
        for s=1:numSampPatients
            tmp2(:,s) = H0_xTest(randi(numPatientSets,numPatientSets,1),s);
        end
        H0_nb_abnormal_patients_ci(:,b) = sum(tmp2<ci(1),1)+sum(tmp2>ci(2),1);
        H0_nb_abnormal_patients_ml(:,b) = sum(tmp2<ci_mlb(1),1)+sum(tmp2>ci_mlb(2),1);
    end
    
    % now we have (for each patient sample tested) a distribution of the
    % classification knowing that there is no difference (since we use H0
    % data)
    for s=1:numSampPatients
        pb = sum(nb_abnormal_patients_ci{i}(s)>H0_nb_abnormal_patients_ci(s,:)) / 600;
        p_value = min(pb,1-pb);
        h_ci(s) = p_value <= alpha;
        
        pb = sum(nb_abnormal_patients_ml{i}(s)>H0_nb_abnormal_patients_ml(s,:)) / 600;
        p_value = min(pb,1-pb);
        h_ml(s) = p_value <= alpha;
    end
    resultsClassif_ci(i) = 100*sum(h_ci)/numSampPatients;
    resultsClassif_ml(i) = 100*sum(h_ml)/numSampPatients;
    
    
    %% -- t-test --
    
    % for the ttest you have to compare with the controls and not with mu
    % this a 2 samples t-tests (the denominator accounts for different
    % variances)
    %     [h, p, ci_ttest] = ttest(xTest', mu, 'Alpha', alpha, 'Tail', 'both'); % h==1 when H0 is rejected (i.e. H0: mean==mu)
    %     resultsTtest(i) = 100*sum(h)/length(h);
    
    [h, p, ci_ttest] = ttest2(repmat(xControl',1,numSampPatients),xTest, ...
        'Alpha', alpha, 'Tail', 'both','vartype','unequal'); % h==1 when H0 is rejected (i.e. H0: mean==mu)
    resultsTtest(i) = 100*sum(h)/numSampPatients;
    
end
%% -- Plot rejection percentage in terms of Delta --
figure;
plot(delta, resultsClassif_ci, 'r','LineWidth',3); 
hold on; 
plot(delta, resultsClassif_ml, 'b','LineWidth',3); 
plot(delta, resultsTtest, 'k--')
legend('Classification ci', 'Claasification ml','t test');
title(['% of rejection in function of Delta in Classification and t-test (alpha = ' num2str(alpha) ')'],'FontSize',16);
xlabel('Delta');
ylabel('% of rejection');
grid on


