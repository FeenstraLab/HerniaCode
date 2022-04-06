
/services/tools/R-3.2.1/bin/Rscript --vanilla /home/projects/ssi_gen1/data/UKbiobank/genotype/scripts/step1_fitNULLGLMM.R \
	--plinkFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/step1_SAIGE.ch1-22_clean \
	--phenoFile=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/hernias_pheno_plus_covariates_v_flashPCA2_fanta4_newINGbroad.txt \
	--phenoCol=newINGbroad \
	--covarColList=sex,yob,PC1,PC2,PC3,PC4,PC5,PC6 \
	--sampleIDColinphenoFile=IID \
	--traitType=binary \
	--outputPrefix=/home/projects/ssi_gen1/data/UKbiobank/genotype/SAIGE/newINGbroad_flashPCA2 \
	--nThreads=22 \
	--LOCO=TRUE

