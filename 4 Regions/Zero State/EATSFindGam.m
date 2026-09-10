% Create a function to find EATS %
% VSI Two-sided 4 Regions Run Sum Xbar Chart for Gamma %
% Zero State Case %

function [EATS] = EATSFindGam(n, K, deltamin, deltamax, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale)

h = 9 ; 

[x, w] = lgwt_table(h) ; 

EATSsum = 0 ; 

for i = 1:h 
    xi = x(i) ; 
    wi = w(i) ; 
    d = ((deltamax - deltamin)/2) * xi + ((deltamax + deltamin)/2) ; 
    ATS = ATSFindGam(n, K, d, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ; 
    EATS = (1/(deltamax - deltamin)) * ATS * wi ; 
    EATSsum = EATSsum + EATS ; 
end

EATS = ((deltamax - deltamin)/2) * EATSsum ; 

end