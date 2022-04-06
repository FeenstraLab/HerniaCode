
/services/tools/R-3.2.1/bin/R --no-save<<EOF
library(MultiPhen)
pheno_cov <- read.table("pheno_and_covariates_unrelated.txt", sep="\t", header=T)
ING <- read.table("ukb_imp_all_chr_unrelated.raw", sep=" ", header=T)
ING_snps <- as.matrix(ING[,7:53])
row.names(ING_snps) <- ING[,1]
ING_pheno <- as.matrix(pheno_cov[,c(2,3,5:14,16)])
row.names(ING_pheno) <- ING[,1]
opts = mPhen.options(c("regression","pheno.input"))
opts$mPhen.variable.selection=TRUE
UMB =mPhen(ING_snps, ING_pheno, phenotypes = all, covariates = c('sex','yob','PC1','PC2','PC3','PC4','PC5','PC6'), opts=opts)
UMB_selection <- data.frame(cbind(rownames(UMB$Res[,,,1]), UMB$Res[,,,1],UMB$Res[,,,2],UMB$Res[,,,3],UMB$maf))
colnames(UMB_selection) <- c("SNPID_EA","DIA_beta","FEM_beta","ING_beta","UMB_beta","VENTR_beta","JointModel_beta","DIA_p","FEM_p","ING_p","UMB_p","VENTR_p","JointModel_p","DIA_N","FEM_N","ING_N","UMB_N","VENTR_N","JointModel_N","MAF")
# Save file 
write.table(UMB_selection, file="ING_PhenoNAs_MultiPhen.txt", row.names = F, sep="\t", quote = F)


