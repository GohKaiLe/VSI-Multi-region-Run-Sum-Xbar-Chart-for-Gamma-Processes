% ATS Optimization for the Two-sided 7 Regions VSI Run Sum Xbar Chart %
% Subject to ATS0 and ASI0 %
% Gamma Distribution %

tic ; 

clearvars; 
clc ; 

% Gamma Distribution Parameters
a_shape = 4 ; % shape parameter 
b_scale = 1 ; % scale parameter

ATS0_target = 370.4 ;
ASI0_target = 1.0 ; 

n = 5 ; % sample size 
delta = 0.50 ; % shift size 

hS = 0.01 ; % short sampling interval 
G = 3 ; % value set by user 

Opt_S1 = 0 ; Opt_S2 = 0 ; Opt_S3 = 0 ; Opt_S4 = 0 ; Opt_S5 = 0 ; Opt_S6 = 0 ; Opt_S7 = 0 ; 
Opt_K = 0 ; Opt_hL = 0 ;

ATS1min = 10000 ; 

out = [] ; 

for S1 = 0:1
    for S2 = S1:3
        for S3 = S2:4 
            for S4 = S3:6
                for S5 = S4:5
                    for S6 = S5: 8
                        for S7 = S6:10
                            [hL, K] = FindhLKVsiRSXbarGam(ATS0_target, ASI0_target, n, hS, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale); 
                            ATS0 = ATSFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ; 
                            ASI0 = ASIFindGam(n, K, 0, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ;
                            if abs(ATS0 - ATS0_target) < 1 
                                ATS1 = ATSFindGam(n, K, delta, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ; 
                                SDTS1 = SDTSFindGam(n, K, delta, hS, hL, G, S1, S2, S3, S4, S5, S6, S7, a_shape, b_scale) ; 
                                out = [out; hS, hL, S1, S2, S3, S4, S5, S6, S7, K, ATS0, ASI0, ATS1, SDTS1] 
                                if ATS1 < ATS1min
                                SDTS1min = SDTS1 ; 
                                ATS1min = ATS1 ;
                                ATS0min = ATS0 ; 
                                ASI0min = ASI0 ; 
                                Opt_hL = hL ; 
                                Opt_K = K ; 
                                Opt_S1 = S1 ; 
                                Opt_S2 = S2 ; 
                                Opt_S3 = S3 ; 
                                Opt_S4 = S4 ; 
                                Opt_S5 = S5 ; 
                                Opt_S6 = S6 ; 
                                Opt_S7 = S7 ; 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

disp(['For Gamma distribution    : Shape (a)  = ', num2str(a_shape,'%0.0f'), ', Scale (b) = ', num2str(b_scale,'%0.0f')])
disp(['Sampling Intervals        : Short (hS) = ', num2str(hS,'%0.2f'), ', Long (hL) = ', num2str(Opt_hL,'%0.2f')])
disp(['The contorl rate (G) is   : ', num2str(G,'%0.0f')])
disp(['The sample size is        : ', num2str(n,'%0.0f')])
disp(['The shift size (delta) is : ', num2str(delta,'%0.2f')])
disp(['The scores are            : Score 1 = ', num2str(Opt_S1,'%0.0f'), ', Score 2 = ', num2str(Opt_S2,'%0.0f'), ', Score 3 = ', num2str(Opt_S3,'%0.0f'),', Score 4 = ', num2str(Opt_S4,'%0.0f')])
disp(['                            Score 5 = ', num2str(Opt_S5,'%0.0f'), ', Score 6 = ', num2str(Opt_S6,'%0.0f'), ', Score 7 = ', num2str(Opt_S7,'%0.0f')])
disp(['The K is                  : ', num2str(Opt_K,'%0.4f')])
disp(['The ATS1 is               : ', num2str(ATS1min,'%0.2f')])
disp(['The SDTS1 is              : ', num2str(SDTS1min,'%0.2f')])
disp(['The ATS0 is               : ', num2str(ATS0min,'%0.2f')])
disp(['The ASI0 is               : ', num2str(ASI0min,'%0.2f')])

toc ; 