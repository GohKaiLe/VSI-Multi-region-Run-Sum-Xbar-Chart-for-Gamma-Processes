% Optimal Parameters for VSI Two-sided 4 Regions RS for Gamma %
% Steady State % 
% Subject to EAATS0 = AATS0 = tau %

tic ;

clearvars; 
clc ; 

% Gamma Distribution Parameters
a_shape = 1 ; % shape parameter 
b_scale = 1 ; % scale parameter

n = 5 ; % sample size 
deltamin = 0.03 ; % minimum shift size 
deltamax = 2.00 ; % maximum shift size 

hS = 0.01 ; % short sampling interval 
G = 3 ; % value set by user 

EAATS0_target = 370.4 ;
ASI0_target = 1.0 ; 

Opt_S1 = 0 ; Opt_S2 = 0 ; Opt_S3 = 0 ; Opt_S4 = 0 ; Opt_K = 0 ; Opt_hL = 0 ;

EAATSmin = 10000 ; 

out = [] ;

for S1 = 0:1
    for S2 = S1:3
        for S3 = S2:5 
            for S4 = S3:10
                [hL, K] = EAATSFindhLKVsiRSXbarGam(EAATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, a_shape, b_scale); 
                EAATS0 = AATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ; 
                ASI0 = ASIFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ;
                if abs(EAATS0 - EAATS0_target) < 1 
                    EAATS = EAATSFindGam(n, K, deltamin, deltamax, hS, hL, G, S1, S2, S3, S4, a_shape, b_scale) ; 
                    out = [out; hS, hL, S1, S2, S3, S4, K, EAATS0, ASI0, EAATS] 
                    if EAATS < EAATSmin
                        EAATSmin = EAATS ;
                        EAATS0min = EAATS0 ; 
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
Final_results = [hS, Opt_hL, Opt_S1, Opt_S2, Opt_S3, Opt_S4, Opt_K, EAATS0min, ASI0min, EAATSmin] 

disp(['For Gamma distribution      : Shape (a)  = ', num2str(a_shape,'%0.0f'), ', Scale (b) = ', num2str(b_scale,'%0.0f')])
disp(['Sampling Intervals          : Short (hS) = ', num2str(hS,'%0.4f'), ', Long (hL) = ', num2str(Opt_hL,'%0.4f')])
disp(['The contorl rate (G) is     : ', num2str(G,'%0.0f')])
disp(['The sample size is          : ', num2str(n,'%0.0f')])
disp(['The range of shift sizes is : Minimum shift = ', num2str(deltamin,'%0.2f'), ', Maximum shift = ', num2str(deltamax,'%0.2f')])
disp(['The scores are              : Score 1 = ', num2str(Opt_S1,'%0.0f'), ', Score 2 = ', num2str(Opt_S2,'%0.0f'), ', Score 3 = ', num2str(Opt_S3,'%0.0f'),', Score 4 = ', num2str(Opt_S4,'%0.0f')])
disp(['The K is                    : ', num2str(Opt_K,'%0.4f')])
disp(['The EAATS is                : ', num2str(EAATSmin,'%0.2f')])
disp(['The EAATS0 is               : ', num2str(EAATS0min,'%0.2f')])
disp(['The ASI0 is                 : ', num2str(ASI0min,'%0.2f')])

toc ; 
