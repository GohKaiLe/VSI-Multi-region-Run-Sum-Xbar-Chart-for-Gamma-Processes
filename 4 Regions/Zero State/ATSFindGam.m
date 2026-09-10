% Create a fucntion to find ATS (Gamma Case) %

function [ATS] = ATSFindGam(n, K, delta, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale)

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

ATS = s' * (inv(I - Q)) * h - s' * h  ;

end