% Optimal Parameters for VSI Two-sided 4 Regions RS for Gamma %
% Zero State % 
% Subject to EATS0 = ATS0 = tau %

tic ;

clearvars; 
clc ; 

% Gamma Distribution Parameters
a_shape = 1.4148 ; % shape parameter 
b_scale = 36.4364 ; % scale parameter

n = 5 ; % sample size 
deltamin = 0.03 ; % minimum shift size 
deltamax = 2.00 ; % maximum shift size 

hS = 0.01 ; % short sampling interval 
G = 3 ; % value set by user 

EATS0_target = 370.4 ;
ASI0_target = 1.0 ; 

Opt_S1 = 0 ; Opt_S2 = 0 ; Opt_S3 = 0 ; Opt_S4 = 0 ; Opt_K = 0 ; Opt_hL = 0 ;

EATSmin = 10000 ; 

out = [] ;

for S1 = 0:1
    for S2 = S1:3
        for S3 = S2:5 
            for S4 = S3:10
                [hL, K] = EATSFindhLKVsiRSXbarGam(EATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, a_shape, b_scale); 
                EATS0 = ATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ; 
                ASI0 = ASIFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ;
                if abs(EATS0 - EATS0_target) < 1 
                    EATS = EATSFindGam(n, K, deltamin, deltamax, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ;
                    out = [out; hS, hL, S1, S2, S3, S4, K, EATS0, ASI0, EATS] 
                    if EATS < EATSmin
                        EATSmin = EATS ;
                        EATS0min = EATS0 ; 
                        ASI0min = ASI0 ; 
                        Opt_hL = hL ; 
                        Opt_K = K ; 
                        Opt_S1 = S1 ; 
                        Opt_S2 = S2 ; 
                        Opt_S3 = S3 ; 
                        Opt_S4 = S4 ; 
                    end
                end
            end
        end
    end
end

disp('The optimal charting parameters are:')
Final_results = [hS, Opt_hL, Opt_S1, Opt_S2, Opt_S3, Opt_S4, Opt_K, EATS0min, ASI0min, EATSmin] 

disp(['For Gamma distribution      : Shape (a)  = ', num2str(a_shape,'%0.0f'), ', Scale (b) = ', num2str(b_scale,'%0.0f')])
disp(['Sampling Intervals          : Short (hS) = ', num2str(hS,'%0.4f'), ', Long (hL) = ', num2str(Opt_hL,'%0.4f')])
disp(['The contorl rate (G) is     : ', num2str(G,'%0.0f')])
disp(['The sample size is          : ', num2str(n,'%0.0f')])
disp(['The range of shift sizes is : Minimum shift = ', num2str(deltamin,'%0.2f'), ', Maximum shift = ', num2str(deltamax,'%0.2f')])
disp(['The scores are              : Score 1 = ', num2str(Opt_S1,'%0.0f'), ', Score 2 = ', num2str(Opt_S2,'%0.0f'), ', Score 3 = ', num2str(Opt_S3,'%0.0f'),', Score 4 = ', num2str(Opt_S4,'%0.0f')])
disp(['The K is                    : ', num2str(Opt_K,'%0.4f')])
disp(['The EATS is                 : ', num2str(EATSmin,'%0.2f')])
disp(['The EATS0 is                : ', num2str(EATS0min,'%0.2f')])
disp(['The ASI0 is                 : ', num2str(ASI0min,'%0.2f')])

toc ; 
