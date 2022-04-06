# QQ plot of inguinal hernia, based on a file with the p-values sorted largest to smallest 

z <-read.table('INGallppadjsortedREV.txt',head=T)
N <-length(z$padj)
z$Pp <- c(N:1)/N-1/(2*N)
png("qqINGallnottruncP.png")
plot(-log10(z$Pp),-log10(z$Pval),xlab="-log10 (expected P)", ylab="-log10 (observed P)",pch=20,cex=0.5,col="red",main="Inguinal")
points(-log10(z$Pp),-log10(z$padj),pch=20,cex=0.5,col="blue")
abline(0,1,col="black")dev.off()

q()

