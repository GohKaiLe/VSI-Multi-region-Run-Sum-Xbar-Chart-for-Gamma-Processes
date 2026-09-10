
function ceq = Constraints(x, ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, a_shape, b_scale)

hL = x(1) ;
K = x(2) ;

ATS0 = ATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ; 
ASI0 = ASIFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ; 

% Equality constraints
ceq(1) = ATS0 - ATS0_target;
ceq(2) = ASI0 - ASI0_target;

end