/services/tools/R-3.2.1/bin/Rscript --vanilla /home/projects/ssi_gen1/data/UKbiobank/genotype/scripts/step2_SPAtests.R \
	--bgenFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/EGAD00010001474/ukb_imp_chr1_v3.bgen \
	--bgenFileIndex=/home/projects/ssi_gen1/data/UKbiobank/genotype/EGAD00010001474/ukb_imp_chr1_v3.bgen.bgi \
    --minMAC=20 \
    --sampleFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/samplefileforbgen.txt \
    --GMMATmodelFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/newINGbroad_flashPCA2.rda \
    --varianceRatioFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/newINGbroad_flashPCA2.varianceRatio.txt \
    --SAIGEOutputFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/chr1.newINGbroad_flashPCA2.bgen.txt \
    --numLinesOutput=10000 \
    --IsOutputAFinCaseCtrl=TRUE

