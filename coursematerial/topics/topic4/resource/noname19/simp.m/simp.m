function q = simp(y,h)

% Parameters
%   y = signal to be integrated over
%   h = step-size

m = length(y)-1; % m = # of sub-intervals (must be even), m+1

if (m/2)~=floor(m/2)
    %disp('Careful...m must be even');
    m = m-1;
end

v = 2*ones(m+1,1);
v2=2*ones(m/2,1);

% length(y)
% size(v)
% size(v2)
% pause;

v(2:2:m) = v(2:2:m)+v2;
v(1) = 1; v(m+1) = 1;

% size(v)
% size(y(1:m+1))
% pause;

q = y(1:m+1)*v;
q=q*h/3;

% If the # of pairs of sub-intervals is not even, use trapezoidal rule
% for last interval

if m+1 < length(y) 
    q = q + (y(end)+y(end-1))*h/2;
end

%q