% About: STATISTICS FOR NEUROSCIENCE RESEARCH

% 9.073 / HST 460

% Class 5: Bayesian Methods 

% Instructor: Dr. Emery N. Brown

%
%-------------------------------------------------------------------------% 
% Code Summary: Rejection Sampling
%
% We assume that a posterior distribution (=3*theta^2) is prescribed.
% We draw samples from this distribution using Rejection Sampling Method. 
% We also determine moments using these samples.
%-------------------------------------------------------------------------% 
%
% Coded by: Sourish Chakravarty (sourish.chakravarty@gmail.com)
% Last updated: 02/27/2017
%
clc 
clear all
close all

% Consider a posterior pdf given by pdf(theta) = 3*theta^2 when theta takes
% values within [0,1] and the function is zero elsewhere.
vec_theta = linspace(0, 1, 100);
pdf_post = 3*vec_theta.^2;% You may verify it has the properties of a pdf


% Now we want to use Rejection Sampling method to generate samples from
% this distribution and also determine, say, its
% first moment, second moment and expectation of a complicated function,
% say, sin(theta*pi). Towards this goal, we need a proposal distribution 
% whose support encompasses the posterior pdf and also there exists a
% constant C such that for all values of theta, 
%(pdf( theta | x ) /h(theta)/C <=1. For the current example, we can use a
%uniform distribution and we can also verify that for C = 3, the required
%condition is satisfied.
C = 3;
fg(1) = figure; hold on
plot( vec_theta, pdf_post, 'b' , 'LineWidth', 2);
plot( vec_theta, ones(size(vec_theta)), 'r' , 'LineWidth', 2);
plot( vec_theta, C*ones(size(vec_theta)), 'k--', 'LineWidth', 2 );
ylim([0,6])
xlabel('theta'), ylabel('pdf')
str_legend = {'Posterior p.d.f.: p(theta|x) = 3*theta^2',...
    'Proposal p.d.f.: h(theta)',...
    'C*h(theta)'};
legend(str_legend, 'Location', 'NorthWest')

% 
N_samp = 1000;
num = [0;0;0];
jcnt = 0;
ncnt = 0;
while 1
    ncnt = ncnt + 1;
    % 1) Draw from Proposal Distribution (A Uniform distribution for this
    % example)
    lambda = rand(1);
    
    % 2) Draw a number from a Uniform Distribution
    u = rand(1);
    
    h_theta = 1;% Calculate the proposal pdf value
    pdf_post = 3*lambda^2;% Calculate the posterior pdf
    % 3) Check condition whether to accept or reject candidate sample
    if u <= pdf_post/(h_theta*C) % Then Accept
        jcnt = jcnt + 1; 
        % Therefore lambda can be considered as a sample drawn from the
        % posterior pdf.
        % 4) Compute g(theta)
        g_theta =[lambda; lambda^2; sin(lambda*pi)];
        % 5) Add contributions to numerator and denominator
        num = num + g_theta;
        if jcnt == N_samp
            break
        end
    end
end
% 6) Calculate Estimates: 
g_theta_est = num/jcnt %[Estimated First Moment; Estimated Second Moment; 
% Expectation of sin(theta*pi)] 
EfficiencyOfAlgorithm = jcnt/ncnt
%--------------------- End of Code ---------------------------------------%