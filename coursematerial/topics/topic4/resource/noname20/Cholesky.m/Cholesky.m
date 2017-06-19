function L = Cholesky(v, p) 
% Wasim Malik 
% Cholesky.m 
% Compute the Cholesky factor L (lower triangular matrix) for the AR(p) 
% process given in vector v 

	N 				= length(v); 									% No of samples of the AR process 
	Linv 			= diag(ones(1, N)); 					% The matrix L^-1 
	Dinv_d 		= zeros(1, N); 								% The diagonal elements of matrix D^-1 
	Dinv_d(1) = 1/xcov(v, 0, 'unbiased');  	% First entry r_xx[0] 
	
	% A matrix: the upper left corner block of Linv 
	for r = 1 : p
		rind = r + 1; 										% Row index for Linv matrix to write on 
		[coef, nvar] = arburg(v, r); 		 	% Estimate AR 
		coef = -coef(2:end); 						 	% Convert z-domain coefs to time-domain coefs of MA(inf) form 
		Linv(rind, 1:r) = -fliplr(coef); 	% Populate first (ARord+1) rows of Linv 
		Dinv_d(rind) = 1/nvar; 						% and first (ARord+1) entries of Dinv diagonal vector 
	end 
	
	% The rest of Linv 
	for r = 1 : (N - p - 1)
		rind = r + p + 1; 											% Row index for Linv matrix to write on 
		Linv(rind, r+1 : r+p) = -fliplr(coef); 	% Write the AR coefficients for order ARord 
	end 
	
	Dinv_d(p+2 : end) = Dinv_d(p+1); 
	L = (Linv' * diag(sqrt(Dinv_d)))'; 

end % function L = Cholesky(v, p)