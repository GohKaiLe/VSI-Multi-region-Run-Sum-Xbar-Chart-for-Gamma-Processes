% Create a fucntion to find the probability matrix Q (Gamma Case) %

function [Q, R1, u] = QMatrixFindGam(n, K, delta, S1, S2, S3, S4, a_shape, b_scale)

a1 = S1 ; % Score region 1 (above mu)
a2 = S2 ; % Score region 2 (above mu)
a3 = S3 ; % Score region 3 (above mu)
a4 = S4 ; % Score region 4 (above mu)

b1 = -a1 ; % Score region 1 (below mu)
b2 = -a2 ; % Score region 2 (below mu)
b3 = -a3 ; % Score region 3 (below mu)
b4 = -a4 ; % Score region 4 (below mu)

% Find transient-probability matrix Q 

B = zeros(4, 1) ;
B(1, 1) = b4 ;
B(2, 1) = b3 ;
B(3, 1) = b2 ;
B(4, 1) = b1 ;

A = zeros(4, 1) ;
A(1, 1) = a1 ;
A(2, 1) = a2 ;
A(3, 1) = a3 ;
A(4, 1) = a4 ;

% Probability for each region 

mu = a_shape * b_scale ; % mean for x bar
var = a_shape * b_scale ^ 2 ; % variance for x bar

UCL1 = mu + K * sqrt (var/n) - delta * sqrt(var) ;
UCL2 = mu + 2 * K * sqrt (var/n) - delta * sqrt(var) ;
UCL3 = mu + 3 * K * sqrt (var/n) - delta * sqrt(var) ;
CL = mu - delta * sqrt(var) ;
LCL1 = mu - K * sqrt (var/n) - delta * sqrt(var) ;
LCL2 = mu - 2 * K * sqrt (var/n) - delta * sqrt(var) ;
LCL3 = mu - 3 * K * sqrt (var/n) - delta * sqrt(var) ;

if LCL1 < 0 
    LCL1 = 0 ;
end

if LCL2 < 0 
    LCL2 = 0 ;
end

if LCL3 < 0 
    LCL3 = 0 ;
end

p3 = gamcdf(UCL3, n * a_shape, b_scale / n) - gamcdf(UCL2, n * a_shape, b_scale / n) ;
p2 = gamcdf(UCL2, n * a_shape, b_scale / n) - gamcdf(UCL1, n * a_shape, b_scale / n) ;
p1 = gamcdf(UCL1, n * a_shape, b_scale / n) - gamcdf(CL, n * a_shape, b_scale / n) ;
p_1 = gamcdf(CL, n * a_shape, b_scale / n) - gamcdf(LCL1 , n * a_shape, b_scale / n) ; 
p_2 = gamcdf(LCL1 , n * a_shape, b_scale / n) - gamcdf(LCL2 , n * a_shape, b_scale / n) ; 
p_3 = gamcdf(LCL2 , n * a_shape, b_scale / n) - gamcdf(LCL3 , n * a_shape, b_scale / n) ;

rows = 100 ;
R = zeros(rows, 18) ;
S = zeros(rows, 2) ;
S(1, 1) = 0 ;
S(1, 2) = 0 ;
u = 1; % Start from column 1

for a = 1:rows
    R(a, 1) = S(a, 1) ;
    R(a, 2) = S(a, 2) ;
    i = 3 ;
    for j = 3:10
        if mod(j, 2) == 1
            R(a, j) = 0 ;
        else
            R(a, j) = R(a, 2) + B(j-i, 1) ;
            i = i + 1 ;
        end
    end
    i = 10 ;
    for j = 11:18
        if mod(j, 2) == 0
            R(a, j) = 0 ;
        else
            R(a, j) = R(a, 1) + A(j-i, 1) ;
            i = i + 1 ;
        end
    end
    for j = 3:10
        if mod(j, 2) == 0
            if R(a, j) > b4
                ww = 0 ;
                for w = 1:rows
                    if R(a, j) ~= S(w, 2)
                        ww = ww + 1 ;
                    end
                end
                      if ww == rows
                        S(u+1, 1) = 0 ;
                        S(u+1, 2) = R(a,j) ;
                        u = u + 1 ;
                        ww = 0;
                      end
            end
        end
    end
    for j = 11:18
        if mod(j, 2) == 1
            if R(a, j) < a4
            ww = 0 ;
            for w = 1:rows
                if R(a, j) ~= S(w, 1)
                    ww = ww + 1 ;
                end
            end
                  if ww == rows
                    S(u+1, 1) = R(a, j) ;
                    S(u+1, 2) = 0 ;
                    u = u + 1 ;
                    ww = 0 ;
                  end
            end
        end
    end
end

  M = zeros(u, rows) ;

  for i = 1:u
    M(i, i) = 1;
  end

  R1 = M * R ;
  temp1 = zeros(1, 18) ;

  for t = 2:(u-1)
    for v = (t+1):u
        if R1(t, 2) > R1(v, 2)
            temp = R1(t, 2) ;
            R1(t, 2) = R1(v, 2) ;
            R1(v, 2) = temp ;
            for i = 1:18
                if i ~= 2
                    temp1(1, i) = R1(t, i) ;
                    R1(t, i) = R1(v, i) ;
                    R1(v, i) = temp1(1, i) ;
                end
            end
        end
    end
  end
  for t = 2:(u-1)
    for v = (t+1):u
        if R1(t, 1) > R1(v, 1)
            temp = R1(t, 1) ;
            R1(t, 1) = R1(v, 1) ;
            R1(v, 1) = temp ;
            for i = 1:18
                if i ~= 1
                    temp1(1, i) = R1(t, i) ;
                    R1(t, i) = R1(v, i) ;
                    R1(v, i) = temp1(1, i) ;
                end
            end
        end
    end
  end

  R2 = zeros(u, 9) ;

  for b = 1:u
    c = 1 ;
    for a = 3:2:18
        d = 0 ;
        for row = 1:u
            if R1(b, a) == R1(row, 1) && R1(b, a + 1) == R1(row, 2)
                R2(b, a - c) = row ;
                d = 1 ;
            end
        end
            if d ~= 1
                R2(b, a-c) = u+1;
            end
                c = c + 1 ;
    end
  end
  for a = 1:u
    R2(a, 1) = a;
  end

  Q = zeros(u, u) ;
  Pr = zeros(1, 6) ;

  Pr(1, 1) = p_3 ;
  Pr(1, 2) = p_2 ;
  Pr(1, 3) = p_1 ;
  Pr(1, 4) = p1 ;
  Pr(1, 5) = p2 ;
  Pr(1, 6) = p3 ;

  for a = 1:u
    for f = 1:u
        prob = 0;
        for b = 3:8
            if R2(a, b) == f
                prob = prob + Pr(1, b - 2) ;
            end
        end
              Q(a, f) = prob;
    end
  end

end