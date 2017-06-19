% About: STATISTICS FOR NEUROSCIENCE RESEARCH

% 9.073 / HST 460

% Class 5: Bayesian Methods 

% Instructor: Dr. Emery N. Brown

%
%-------------------------------------------------------------------------% 
% Code Summary: Gaussian Approximation & Importance Sampling
%
% Part (A) We use N number of I.I.D. values drawn from a Poisson
% distribution as our observation. Then using a Gamma distribution as prior
% we generate a Gamma posterior pdf on the unknown distribution parameter
% (of the distribution from which the samples are drawn) that we are trying 
% to estimate. In the matlab code below we plot the Observations, 
% the Prior pdf, the Likelihood and the exact Posterior. We finally plot a 
% Gaussian approximation of the exact Posterior pdf. You may change the
% value of N to view how it affects the relevant pdfs and estimates.
%
% Part (B) We want to determine the first and second moments of the
% exact posterior distribution via Importance Sampling using the 
% Gaussian approximation as the Importance Distribution. 
% The use of the approximation of the posterior pdf as the
% Importance pdf is a convenient choice. We could have used any other pdf
% as long as the support of this candidate pdf will match that of the
% posterior distribution from which we want to draw samples. A reasonable
% choice is a pdf which resembles the posterior distribution and is
% easy to draw from. Hence, our current choice is a reasonable one.
%-------------------------------------------------------------------------% 
%
% Coded by: Sourish Chakravarty (sourish.chakravarty@gmail.com)
% Last updated: 02/27/2017
%
clc 
clear all
close all
%% Generating observation with I.I.D. Poisson random variables
N = 15; % Number of observations
lambda = 5;% pdf parameter to draw observations
x_samp = poissrnd(lambda, [N, 1]);%x_samp contains non-negative integers
fg(1) = figure; hold on
[nfreq, xout] = hist(x_samp, [min(x_samp):max(x_samp)] );
plot([xout; xout], [zeros(1, length(xout)); nfreq/sum(nfreq)],...
    'b', 'LineWidth',5 );%Plot the observations
xlabel('Observation, x');
ylabel('pdf');
%-- From here onwards we shall assume that <lambda> used to generate the
%samples is unknown and we only have the observation vector <x_samp>. 

%% Part A: Gaussian approximation of posterior pdf
% Since <lambda> is unknown, we assume that we have some prior knowledge 
% in terms of a probability distribution on <lambda>: the prior
% distribution. Since <lambda> is a non-negative quantity and the sampling
% distribution is Poisson, therefore a suitable conjugate prior for 
% <lambda> is a % Gamma distribution(kAlpha, kBeta) where the 
% kAlpha (shape parameter) and kBeta (scale parameter) are
% hyper-parameters that characterize the Gamma pdf for <lambda>, such that 
% Mean = kAlpha/kBeta, Variance = kAlpha/(kBeta^2)
kBeta = 1;
kAlpha = 10;
vecL = linspace(1, kAlpha*2, 100);% Create a grid in <lambda> axis
pdf_prior = gampdf(vecL, kAlpha, 1/kBeta);% calculate prior
% Note: the definition of the scale parameter in Matlab documentation is
% inverse of the scaling term in the Gamma p.d.f. expression as defined in
% the class notes.


% To calculate the likelihood function the summary statistic we use here 
% is: Sn = sum(x_samp)
Sn = sum(x_samp);
pdf_likelihood = (vecL.^Sn).*exp(-N*vecL)/prod(factorial(x_samp));
% This is the joint pdf of the observations. Since we have substituted the
% observation in this expression, this becomes the likelihood which is 
% a function of the distribution parameter <lambda>.

% Note: since the values in pdf_likelihood are very small numbers, we
% divide everywhere by the maximum value so that its shape can be compared
% with the prior and posterior pdf-s.

% As seen in the notes, the exact posterior distribution on <lambda> 
% conditioned on the given data, 
% for this scenario, is given by another Gamma pdf characterized 
% by a shape parameter = (Sn + kAlpha) and the 
% scale parameter = (kBeta + N). Now we calculate the exact
% posterior based on the above assumptions. 
pdf_post_exact = gampdf(vecL, Sn + kAlpha, 1/(kBeta + N) );


%-- Gaussian approximation
% Here we use the Gaussian approximation where the posterior is linearized
% about the MAP estimate. (See notes for derivation)
mean_gauss = (Sn + kAlpha -1)/(kBeta + N + 1);
var_gauss = (Sn + kAlpha -1)/(kBeta + N + 1)^2;
pdf_post_approx = normpdf(vecL,mean_gauss,var_gauss);

%-- Plot the pdfs
fg(2) = figure; hold on
plot(vecL, pdf_prior, 'r','LineWidth', 2)
plot(vecL, (1/max(pdf_likelihood))*pdf_likelihood, 'b','LineWidth', 2)
plot(vecL, pdf_post_exact, 'm','LineWidth', 2)
plot(vecL, pdf_post_approx, 'k--','LineWidth', 2)
xlabel('Lambda');
ylabel('pdf');
cell_legend = {'Prior p.d.f.',...
    ['Likelihood/(', num2str(max(pdf_likelihood)),')'],...
    'Exact Posterior p.d.f.',...
    'Gaussian Approximation to Posterior p.d.f.'};
legend(cell_legend, 'Location', 'NorthEast');

%% Part B: Importance Sampling 
% We use the Gaussian approximation to the posterior pdf as our
% Importance Distribution to estimate the first and second moments of the
% posterior distribution;
N_samp = 1000;
num = [0;0];
den = 0;
for ncnt = 1:N_samp
    % 1) Draw from Importance Distribution.
    lambda = mean_gauss + sqrt(2)*var_gauss*randn(1);
    % Again note that the lambda has not been drawn from the posterior pdf.
    % Since our goal is to get reliable estimates of integrals 
    % (moments are defined in terms of expectation operation, 
    % hence they are integrals for continuous random variables) we can get
    % by with drawing samples from a closely resembling
    % distribution (the Importance Distribution) and weighting the 
    % integrand values appropriately to eventually arrive at the right 
    % result for the expectation terms.
    %
    % 2) Compute w(theta|x) and g(theta)
    g_theta = [lambda; lambda^2];
    h_theta =  normpdf(lambda,mean_gauss,var_gauss) ;
    pdf_likelihood = (lambda.^Sn).*exp(-N*lambda)/prod(factorial(x_samp));
    pdf_prior = gampdf(lambda, kAlpha, 1/kBeta);
    w_theta = pdf_likelihood*pdf_prior/h_theta ;
    % 3) Add contributions to numerator and denominator
    num = num + g_theta*w_theta;
    den = den + w_theta;
end
% 4) Final Estimate is numerator/denominator
g_theta_est = num/den % [Estimated First Moment; Estimated Second Moment] 
%--------------------- End of Code ---------------------------------------%