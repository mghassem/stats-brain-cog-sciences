

%generate 600 trials of binary data from the
%logistic equation 

p0 = 0.25; pf = 0.95; 

t=[1:600];
psim(t) = p0 + (pf - p0)./(1.0+exp(-0.1*(t-300)));

prand = rand(1, 600);

behavior = prand<psim;

save behavior behavior;
    