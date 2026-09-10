
function ceq = EATSConstraints(x, EATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale)

    hL = x(1);
    K = x(2);

    EATS0 = ATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;
    ASI0 = ASIFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;

    % Check for Inf/NaN
    if isnan(EATS0) || isinf(EATS0) || isnan(ASI0) || isinf(ASI0)
        ceq = [1e10; 1e10];  % Penalize invalid values
    else
        ceq = [EATS0 - EATS0_target; ASI0 - ASI0_target];
    end

end