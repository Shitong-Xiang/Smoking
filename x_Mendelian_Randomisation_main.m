%% Prepare GWAS and PRS
clear;clc

path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/L_vmPFC';
bfile = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/gene_imputation/imputed_plink_qc_merge_rs/imagen_merge_all_unique_filt_rs';
pheno = fullfile(path,'IMAGEN_L_vmPFC.txt');
covar = fullfile(path,'IMAGEN_Cova.txt');
Out = fullfile(path,'GWAS_Out');
unix(['/share/inspurStorage/home1/ISTBI_data/toolbox/plink2 ' ...
    '--bfile ' bfile ...
    ' --pheno ' pheno ...
    ' --covar ' covar ...
    ' --linear hide-covar --maf 0.01' ...
    ' --covar-variance-standardize' ... %--missing-var-code NA
    ' --out ' Out]);

gwas_file  = fullfile(path,'GWAS_Out.Phenotype.glm.linear'); %% Rename
out = fullfile(path,'imagen_merge_all_unique_filt_rs');
unix(['plink1.9 --bfile ',bfile,...
    ' --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump ',gwas_file,...
    ' --clump-snp-field SNP --clump-field P --out ',out]);

clump = [out,'.clumped'];
Valid_snp = fullfile(path,'GWAS.valid.snp');
unix(['awk ''NR!=1{print $3}'' ',clump,' > ',Valid_snp]);
SNP_p = fullfile(path,'GWAS.SNP.pvalue');
unix(['awk ''{print $3,$12}'' ',gwas_file,' > ',SNP_p]);
range_list = fullfile(path,'range_list');
unix(['echo "0.001 0 0.001" > ',range_list]);
unix(['echo "0.05 0 0.05" >> ',range_list]);
unix(['echo "0.1 0 0.1" >> ',range_list]);
unix(['echo "0.2 0 0.2" >> ',range_list]);
unix(['echo "0.3 0 0.3" >> ',range_list]);
unix(['echo "0.4 0 0.4" >> ',range_list]);
unix(['echo "0.5 0 0.5" >> ',range_list]);
PRS_out = fullfile(path,'PRS');
unix(['plink1.9 --bfile ',bfile,' --score ',gwas_file,...
    ' 3 6 9 header --q-score-range ',range_list,' ',SNP_p,' --extract ',Valid_snp,' --out ',PRS_out]);

%% Check the best threshold
clear;clc

Cova = readtable('/path/to/Covariants.csv');
Pheno = readtable('/path/to/Phenotype.csv');
PRS = dir('/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/L_vmPFC/PRS.*.profile');
for i = 1:length(PRS)
    Score = importdata(fullfile(PRS(i).folder,PRS(i).name));
    Score = Score.data;
    [~,ind1,ind2] = intersect(Pheno.SubID,Score(:,2));
    Data(:,i) = Score(ind2,end);
end

[r,p] = partialcorr(Pheno.Phenotype(ind1),Data,table2array(Cova(ind1,2:end)),'row','complete');
[~,thr_ind] = max(r);

%% MR analysis
clear;clc

output = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/L_vmPFC-Smoking';
pheno1_path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/L_vmPFC';
pheno2_path = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/Smoking';

pheno1_gwas_file = fullfile(pheno1_path,'GWAS_Out.Phenotype.glm.linear');
pheno1_clump = fullfile(pheno1_path,'imagen_merge_all_unique_filt_rs.clumped');
pheno1_Valid_SNP = fullfile(output,'Pheno1_Valid_SNP.txt');
unix(['awk ''{print $3,$5}'' ',pheno1_clump,' > ',pheno1_Valid_SNP]);
pheno1_Valid_SNP = readtable(pheno1_Valid_SNP);
pheno1_best_SNP = pheno1_Valid_SNP;

pheno2_gwas_file = fullfile(pheno2_path,'GWAS_Out.Phenotype.glm.linear');
pheno2_clump = fullfile(pheno2_path,'imagen_merge_all_unique_filt_rs.clumped');
pheno2_Valid_SNP = fullfile(output,'Pheno2_Valid_SNP.txt');
unix(['awk ''{print $3,$5}'' ',pheno2_clump,' > ',pheno2_Valid_SNP]);
pheno2_Valid_SNP = readtable(pheno2_Valid_SNP);
pheno2_best_SNP = pheno2_Valid_SNP;

Threshold = [0.001,0.05,0.1,0.2,0.3,0.4,0.5];
bfile = '/share/inspurStorage/home1/ISTBI_data/IMAGEN/gene_imputation/imputed_plink_qc_merge_rs/imagen_merge_all_unique_filt_rs';

pheno1_best_SNP_Thr = 0.05;
pheno1_best_SNP(pheno1_Valid_SNP.P > pheno1_best_SNP_Thr,:) = [];
[~,ind1,ind2] = intersect(pheno1_best_SNP.SNP,pheno2_Valid_SNP.SNP);
for i = 1:length(Threshold)
    tmp = pheno1_best_SNP;
    tmp(ind1(pheno2_Valid_SNP.P(ind2)<Threshold(i)),:) = [];
    tmp(isnan(tmp.P),:) = [];
    writetable(tmp,fullfile(output,'tmp.txt'),'Delimiter',' ');
    pheno1_Select_SNP = fullfile(output,['Pheno1_Select_Pheno2_Thr_',replace(num2str(Threshold(i)),'.',''),'.txt']);
    unix(['awk ''NR!=1{print $1}'' ',fullfile(output,'tmp.txt'),' > ',pheno1_Select_SNP]);
    unix(['rm ',fullfile(output,'tmp.txt')]);
    Pheno1_PRS = fullfile(output,['Pheno1_PRS_Ex_Pheno2_Thr_',replace(num2str(Threshold(i)),'.','')]);
    unix(['plink1.9 --bfile ',bfile,' --score ',pheno1_gwas_file,' 3 6 9 header --extract ',pheno1_Select_SNP,' --out ',Pheno1_PRS]);
end

pheno2_best_SNP_Thr = 0.1;
pheno2_best_SNP(pheno2_Valid_SNP.P > pheno2_best_SNP_Thr,:) = [];
[~,ind1,ind2] = intersect(pheno2_best_SNP.SNP,pheno1_Valid_SNP.SNP);
for i = 1:length(Threshold)
    tmp = pheno2_best_SNP;
    tmp(ind1(pheno1_Valid_SNP.P(ind2)<Threshold(i)),:) = [];
    tmp(isnan(tmp.P),:) = [];
    writetable(tmp,fullfile(output,'tmp.txt'),'Delimiter',' ');
    pheno2_Select_SNP = fullfile(output,['Pheno2_Select_Pheno1_Thr_',replace(num2str(Threshold(i)),'.',''),'.txt']);
    unix(['awk ''NR!=1{print $1}'' ',fullfile(output,'tmp.txt'),' > ',pheno2_Select_SNP]);
    unix(['rm ',fullfile(output,'tmp.txt')]);
    Pheno2_PRS = fullfile(output,['Pheno2_PRS_Ex_Pheno1_Thr_',replace(num2str(Threshold(i)),'.','')]);
    unix(['plink1.9 --bfile ',bfile,' --score ',pheno2_gwas_file,' 3 6 9 header --extract ',pheno2_Select_SNP,' --out ',Pheno2_PRS]);
end

Cova = readtable('/path/to/Covariants.csv');
Beh1 = readtable('/path/to/Phenotype1.txt');
Beh2 = readtable('/path/to/Phenotype2.txt');

[~,ind1,ind2] = intersect(Cova.SubID,Beh2.SubID);
Cova = Cova(ind1,:);
Beh1 = [Beh1.BL_L_vmPFC(ind1),Beh1.FU2_L_vmPFC(ind1)];
Beh2 = [Beh2.BL_Smoking(ind2),Beh2.FU2_Smoking(ind2)];

PRS1 = dir('/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/L_vmPFC-Smoking/Pheno1*.profile');
for i = 1:length(PRS1)
    Score = importdata(fullfile(PRS1(i).folder,PRS1(i).name));
    Score = Score.data;
    [~,ind1,ind2] = intersect(Cova.SubID,Score(:,2));
    Data(:,i) = Score(ind2,end);
end
[r1,p1] = partialcorr(Beh2(ind1,:),Data,table2array(Cova(ind1,2:end)),'row','complete');
Stat1 = [r1(2,:)',p1(2,:)'];

PRS2 = dir('/share/inspurStorage/home1/ISTBI_data/IMAGEN/Smoke/Mendelian_randomization/L_vmPFC-Smoking/Pheno2*.profile');
for i = 1:length(PRS2)
    Score = importdata(fullfile(PRS2(i).folder,PRS2(i).name));
    Score = Score.data;
    [~,ind1,ind2] = intersect(Cova.SubID,Score(:,2));
    Data(:,i) = Score(ind2,end);
end
[r2,p2] = partialcorr(Beh1(ind1,:),Data,table2array(Cova(ind1,2:end)),'row','complete');
Stat2 = [r2(2,:)',p2(2,:)'];
