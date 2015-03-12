% scrExampleGamma

close all;
clear;

%% Parameters

n1 = 10000;
n2 = 10000;

a = 10; % Shape
b = 0.75; % Scale
delta = 10;

%% Generate distributions

% See the distribution in comparison to the normal distribution
x = 0:0.0002:20; 
y = gampdf(x, a, b); 
figure; plot(x, y, 'b-'); ylim([0,0.3]); hold on;
plot(x,normpdf(x,7,2), 'r--'); hold off;

[g, g1, g2] = fGenerateBimodGamma(n1, n2, a, b, delta);

h1 = hist(g1, 100);
figure; bar(h1, 'b'); ylim([0,600]); %hold on;
h2 = hist(g2, 100);
figure; bar(h2, 'r'); ylim([0,600]); %hold off;
h3 = hist(g, 100);
figure; bar(h3, 'g'); ylim([0,600]);
% figure; hist([g1 g2], 100);
