#单基因生存分析
Rscript "/home/liuchenglong/Documents/program/SurvivalAnalysis/SurvivalAnalysisV1.R" \
--intgene 'SPP1,akshdk' \
--disease 'LIHC' \
--input "/home/liuchenglong/Documents/TCGA/4.ENSEMBL2SYMBOL/LIHC_tpm_unstranded.txt" \
--outdir /home/liuchenglong/Documents/program/SurvivalAnalysis/ \
--cutoff 'median' \
--prefix 'test' \
--timetype month \
--level None
#双基因组合生存分析
Rscript "/home/liuchenglong/Documents/program/SurvivalAnalysis/SurvivalAnalysisV1.R" \
--intgene 'NRG2,NDRG1' \
--disease 'LUAD' \
--input "/home/liuchenglong/Documents/TCGA/4.ENSEMBL2SYMBOL/LUAD_tpm_unstranded.txt" \
--outdir /home/liuchenglong/Documents/program/SurvivalAnalysis/ \
--cutoff 'median' \
--prefix 'test' \
--timetype month \
--level 'low,high,high,low'
