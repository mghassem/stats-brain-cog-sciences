function  U = TRT(dNt,h,lambda)

% Parameters
%   dNt = discrete-time spike train
%   lambda = rate of the process
%   h = sampling interval of CT process from which dNt obtained
% 
%    N.B: dNt and lambda are assumed to be ROW vectors
%
% Output
%   U = 1 - exp(rescISI), where rescISI are rescaled interspike 
%       intervals 


tempdNt = find(dNt~=0); % DT locations of spikes
n = length(tempdNt);
rescfac = zeros(size(tempdNt)-1); % Stores area of lambda b/w consecutive spikes

% Time-rescaling analysis

for i=0:length(tempdNt)-1

   if i==0
       if length(lambda(1:tempdNt(i+1))) > 2

           rescfac(i+1) = simp(lambda(1:tempdNt(i+1)),h);
       else
           rescfac(i+1) = h*trapz(lambda(1:tempdNt(i+1)));               
       end
    
   else
       
        if length(lambda(tempdNt(i):tempdNt(i+1))) > 2
            rescfac(i+1) = simp(lambda(tempdNt(i):tempdNt(i+1)),h);
        else
            rescfac(i+1) = h*trapz(lambda(tempdNt(i):tempdNt(i+1)));
        end
        
   end
end

U = 1 - exp(-rescfac);
idx1 = find(U==1);
if length(idx1) > 0
    U(idx1) = 0.999;
end
idx0 = find(U==0);
if length(idx0) > 0
    U(idx0) = 0.001;
end