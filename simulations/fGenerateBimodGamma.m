function [g, g1, g2] = fGenerateBimodGamma(n1, n2, a, b, delta)
% This function generates a bimodal distribution. The groups have n1 and n2
% samples, respectively. The first one follows a gamma distribution with
% parameters shape=10 and scale=0.75 (so, its mean will be 7.5 and its VAR
% is 5.625). The second group is an inverse gamma with the same parameters
% and separated by delta from the first group
% 
% INPUT
%   - n1: Number of elements in the first group 
%   - n2: Number of elements in the second group 
%   - a: Shape parameter of the gamma distribution (recommended: 10)
%   - b: scale parameter of the gamma distribution (recommended 0.75)
%   - delta: Separation between the two groups
%
% OUTPUT
%   g: All generated samples
%   g1: Samples of the first group 
%   g2: Samples of the second group
%

% % See the distribution in comparison to the normal distribution
%     x = 0:0.0002:20; 
%     a1= 10; b1=0.75; y1 = gampdf(x,a1,b1); 
%     plot(x, y1, 'b-'); ylim([0,0.3]); hold on;
%     plot(x,normpdf(x,7,2), 'r--')

%     a = 10;  % Shape
%     b = 0.75; % Scale

    g1 = gamrnd(a, b, [1, n1]);
    
    g2_aux =  gamrnd(a, b, [1, n2]);
    g2 = delta+(min(g2_aux) + max(g2_aux) - g2_aux);
   
    g = [g1 g2];
    
%     h1 = hist(g1, 100);
%     bar(h1, 'b'); hold on;
%     h2 = hist(g2, 100);
%     bar(h2, 'r'); hold off;
%     h3 = hist([g1 g2], 100);
%     figure; bar(h3, 'g');
%     figure; hist([g1 g2], 100);

end
