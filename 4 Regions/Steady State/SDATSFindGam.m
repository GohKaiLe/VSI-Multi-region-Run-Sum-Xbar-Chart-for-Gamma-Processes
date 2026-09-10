% Create a fucntion to find SDATS %
% Gamma Distribution % 
% Steady State %

function [SDATS] = SDATSFindGam(n, K, delta, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale)

[Q, R1, u] = QMatrixFindGam(n, K, delta, S1, S2, S3, S4, a_shape, b_scale) ; 

s = zeros(u, 1) ; 
s(1, 1) = 1 ; 
I = eye(u) ; 

H = S4 ; 
h = ones(u, 1) ;  
temp2 = R1(:, 1) ; % extract all possible upper cumulativen scores
temp3 = R1(:, 2) ; % extract all possible lower cumulativen scores

for i = 1:u
    if (H/G <= temp2(i)) && (temp2(i) < H)
        h(i) = hS ;
    else
        h(i) = hL ;
    end 
    if (-H < temp3(i)) && (temp3(i) <= -H/G)
        h(i) = hS ;
    end
end

one = ones(u, 1) ; 

s0 = ((inv(I - Q')) * s ) / (one' * (inv(I - Q')) * s) ;

alpha = zeros(u, 1) ; 

for i = 1:u 
    alpha(i) = (s0(i) * h(i)) / (s0' * h) ; 
end

E = diag(h) ; 
V1 = (2 * (inv(I - Q))) - I ;  
V2 = (1/2) * E - (1/4) * h * alpha' ;   

SDATS = sqrt(alpha' * (V1 * V2 * V1 - (1/6) * E) * h) ;

end