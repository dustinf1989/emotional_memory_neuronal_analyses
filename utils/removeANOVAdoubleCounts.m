function anovaResX = removeANOVAdoubleCounts(anovaResX)

for i = 1:length(anovaResX) % Would need to eliminate stubjects here
        iY = (anovaResX{i}.p(3,:) < 0.05);
        anovaResX{i}.p(1,iY) = 1;
        anovaResX{i}.p(2,iY) = 1;
    
        anova_p_sh_noX = anovaResX{i}.p_sh;
        iX = (anovaResX{i}.p_sh(:,:,3) < 0.05); 
        anova_p_sh_noX1 = anova_p_sh_noX(:,:,1);
        anova_p_sh_noX2 = anova_p_sh_noX(:,:,2);
        anova_p_sh_noX1(iX) = 1;
        anova_p_sh_noX2(iX) = 1;
        anovaResX{i}.p_sh(:,:,1) = anova_p_sh_noX1;
        anovaResX{i}.p_sh(:,:,2) = anova_p_sh_noX2;
        paths.pics = 'D:\Zurich_IAPS_micro\results\combinato_Encoding_30Jan2023\withoutDoubleCounting\';
        %  'D:\Zurich_IAPS_micro\results\combinato_Encoding_30Jan2023\AllPics_enc_29Nov23_0.25mfr_50trials_posOnly_4anovas_10000sh_manuK_best5Only\'
end
    
