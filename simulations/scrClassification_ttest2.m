% In this script the comparison between a classification and a t-test with
% normally distributed univariate samples will be compared. Concretely,
% Control and Test sets sampled from a Normal distribution with given  mean
% and std (in the Parameters section) is used.
%
% a) Classification: The control set is used to create a bootstrap CI
% of the means or maximim likelihood CI. After that, the classification
% consists in looking if each patient falls into the CI or not. A bootstrap
% under H0 is then performed to compute the significance of the
% classification result.
%
% b) t-test: assess whether the mean of the test sets can be assumed
% to be equal to mu or not.
%
% In this second simulation we add 1 loop to investifgate the dependency on
% the control sample size since this will influcence the CI range.


clear all
clc

%% -- Parameters --

% mu_control = 130;
% mu_pat = 150;
mu = 130;
sigma = 20;
delta = 0:20;
numDeltas = length(delta);

SampControl = [10 50 100:100:1000];
numSampControl = length(SampControl);
numPatientSets = 200;
numSampPatients = 1000;
nBoot = 5000;
alpha = 0.05;

resultsClassif = zeros(numSampControl, numDeltas);
resultsTtest = zeros(numSampControl, numDeltas);

%% -- Generate control and test samples --
% The control and tests are samples from a Normal distribution with mean mu
% (and mu+alpha in the case of the tests) and std sigma. First, the control
% group is generated, and the bootstrap CI is computed. Then, the whole
% process of generating the control sets, and classification / t-test is
% done once per delta, storint the results in the vectors resultsClassif
% and resultsTtest

%Sample control set
Control = normrnd(mu,sigma, [1, 10000]);

for s=1:numSampControl
    xControl = Control(randi(1000,SampControl(s),1)); % pick SampControl(s) controls
    
    out = prctile(xControl,[0:100]);
    ci = [out(3) out(98)];
    
    nb_abnormal_controls(s) = (sum(xControl<ci(1))+sum(xControl>ci(2))) / SampControl(s) * 100;

    for i=1:numDeltas
        fprintf('running Control N=%g Delta = %g\n',SampControl(s),delta(i));
        xTest =  normrnd(mu+delta(i), sigma, [numPatientSets, numSampPatients]);
        
        %% -- Classification by bootstrap CI --
        nb_abnormal_patients{i} = sum(xTest<ci(1),1)+sum(xTest>ci(2),1);
        
        % now we use bootstrap to evaluate if this number of abnormal patients
        % is signifant (the null hypothesis is controls = patients so we simply
        % do many time the same operation but centering the data)
        H0_xControl = xControl-mean(xControl);
        H0_xTest = xTest - repmat(mean(xTest,1),200,1);
        tmp2 = NaN(numPatientSets,numSampPatients);
        H0_nb_abnormal_patients = NaN(numSampPatients,600);
        for b= 1:600
            % compute a new ci by resampling H0 controls
            tmp = H0_xControl(randi(SampControl(s),SampControl(s),1));
            out = prctile(tmp,[0:100]);
            cib = [out(3) out(98)];
            
            % compute a new classification by resampling H0 patients
            for ss=1:numSampPatients
                tmp2(:,ss) = H0_xTest(randi(numPatientSets,numPatientSets,1),ss);
            end
            H0_nb_abnormal_patients(:,b) = sum(tmp2<cib(1),1)+sum(tmp2>cib(2),1);
        end
        
        % now we have (for each patient sample tested) a distribution of the
        % classification knowing that there is no difference (since we use H0
        % data)
        for ss=1:numSampPatients           
            pb = sum(nb_abnormal_patients{i}(ss)>H0_nb_abnormal_patients(ss,:)) / 600;
            p_value = min(pb,1-pb);
            h_ml(ss) = p_value <= alpha;
        end
        resultsClassif(s,i) = 100*sum(h_ml)/numSampPatients;
        
        
        %% -- t-test --
        
        % for the ttest you have to compare with the controls and not with mu
        % this a 2 samples t-tests (the denominator accounts for different
        % variances)
        %     [h, p, ci_ttest] = ttest(xTest', mu, 'Alpha', alpha, 'Tail', 'both'); % h==1 when H0 is rejected (i.e. H0: mean==mu)
        %     resultsTtest(i) = 100*sum(h)/length(h);
        
        [h, p, ci_ttest] = ttest2(repmat(xControl',1,numSampPatients),xTest, ...
            'Alpha', alpha, 'Tail', 'both','vartype','unequal'); % h==1 when H0 is rejected (i.e. H0: mean==mu)
        resultsTtest(s,i) = 100*sum(h)/numSampPatients;
    end
    classification_rate(s,:,:) = reshape(cell2mat(nb_abnormal_patients),1000,size(nb_abnormal_patients,2));
    clear nb_abnormal_patients
end

figure;
subplot(1,2,1)
plot(delta,resultsClassif,'LineWidth',3); 
grid on; axis([-0.1 20.1 0 101])
title(['Classification power using percentiles (alpha = ' num2str(alpha) ')'],'FontSize',16);
xlabel('Delta'); ylabel('% of rejection'); legend({'10 controls','50 controls','100 controls','200 controls','300 controls', ...
    '400 controls','500 controls','600 controls','700 controls','800 controls','900 controls','1000 controls'},'Location','northwest')

subplot(1,2,2)
plot(delta, resultsTtest,'LineWidth',3); 
grid on; axis([-0.1 20.1 0 101])
title(['T-test power (alpha = ' num2str(alpha) ')'],'FontSize',16);
xlabel('Delta'); ylabel('% of rejection'); legend({'10 controls','50 controls','100 controls','200 controls','300 controls', ...
    '400 controls','500 controls','600 controls','700 controls','800 controls','900 controls','1000 controls'},'Location','northwest')

figure;
plot(delta,squeeze(mean(classification_rate,2))./200.*100,'LineWidth',3); 
grid on; axis([-0.1 20.1 0 101])
title(['Classification rate using percentiles (alpha = ' num2str(alpha) ')'],'FontSize',16);
xlabel('Delta'); ylabel('% of rejection'); legend({'10 controls','50 controls','100 controls','200 controls','300 controls', ...
    '400 controls','500 controls','600 controls','700 controls','800 controls','900 controls','1000 controls'},'Location','northwest')
