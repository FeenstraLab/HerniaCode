# Copyright (C) 2014 Genome Research Ltd / EMBL - European Bioinformatics Institute
#
# Author : Valentina Iotchkova <vi@sanger.ac.uk>
# Author : Matthias Geihs <mg18@sanger.ac.uk>
#
# This file is part of GARFIELD - GWAS analysis of regulatory or functional information enrichment with LD correction.
#
# GARFIELD is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.




#### Read input parameters ####
	#options(echo=TRUE) # if you want see commands in output file
	args <- commandArgs(trailingOnly = TRUE)

	input_file <- args[1] # file with fold enrichment and significance values
	num_perm <- as.integer(args[2]) # number of permutations pre-specified when using the tool
	output_prefix <- args[3] # file prefix for the figures, actual filenames will be "output_prefix"_category, where category is specified in the link_file (e.g. dhs_hotspots etc...)
	plot_title <- args[4] # title for the plot, no title denoted by ""
	filter <-as.integer(args[5]) # filter on the minimum number of variants per treshold
	tr <- -log10(as.numeric(args[6])) ## significance threshold, if 0 then the default threshold below is set up
	if (tr==Inf) tr=-log10(0.05/498)

	rm(args)
	print(paste("Input file: ",input_file,sep=""))
	print(paste("Number of permutations: ",num_perm,sep=""))
	print(paste("Figure prefix: ",output_prefix,sep=""))
	print(paste("Figure title: ",plot_title,sep=""))
	print(paste("Minimum number of variants filter: ",filter,sep=""))
	print(paste("Significance threshold : -log10 Pvalue=",tr,sep=""))

	
#### Read files
	source("garfield-plot-function.R") # load R function to create plot
	input <- read.table(input_file, header=T) # file format: annotation threshold fold_enrichment pvalue ...
	input=input[which(input[,7]>=filter),] # filter variants
	input[which(input[,4]==0),4] = 1/num_perm ## use num_perm to transform 0 pvalues to a finite limit of the -log10 scale
	input[which(input[,4]==(-1)),4] = NA # set -1 pvalues as NAs 
#### Create plots 

	input$Category=as.character(input$Category)
	if (length(which(is.na(input$Category)))>1){ 
		input$Category[which(is.na(input$Category))] = "Custom" 
	}
	input$Category=as.factor(input$Category)

	for (category in as.character(unique(input$Category))) {

		ids = which(input$Category==category)
		thresholds = sort(unique(input[ids,2]))
		thresholdsP = sort(unique(input[ids,2][which(!is.na(input[ids,4]))]))
		annotations = unique(as.character(input[ids,1]))
		if (category %in% c("Genic","Histone_Modifications","Chromatin_States")){
			tissues = as.character(input$Type[ids][match(annotations,input[ids,1])])
			nms = as.character(input$Celltype[ids][match(annotations,input[ids,1])])
			compact=FALSE
			tissue_label = "Feature"
			if (category %in% c("Genic")){
				nms=tissues
			}
		} else if (category %in% c("TFBS","FAIRE","Hotspots","Peaks","Footprints")){
			nms = as.character(input$Celltype[ids][match(annotations,input[ids,1])])
			tissues = as.character(input$Tissue[ids][match(annotations,input[ids,1])])
			compact=TRUE
			tissue_label = "Tissue"
			if (category %in% c("Hotspots","Peaks","Footprints")){
				nms=tissues
			}		
		} else {

			nms = as.character(input$Celltype[ids][match(annotations,input[ids,1])])
			tissues = as.character(input$Tissue[ids][match(annotations,input[ids,1])])
		}


		DATA = matrix(NA, nr=length(thresholds)+1, nc=length(annotations))
		DATA_p = matrix(NA, nr=length(thresholdsP)+1, nc=length(annotations))
		for (j in 1:length(annotations)){
			for (i in 1:length(thresholds)){
				ii = which(input[,1]==annotations[j] & input[,2]==thresholds[i])
				DATA[i,j] = input[ii,3]
			}
			for (i in 1:length(thresholdsP)){
				ii = which(input[,1]==annotations[j] & input[,2]==thresholdsP[i])
				DATA_p[i,j] = input[ii,4]
			}
		}
		DATA[length(thresholds)+1,] = DATA_p[length(thresholdsP)+1,] = 1
		DATA_p = -log10(DATA_p)

		ann.col = colorRampPalette(c("tomato","skyblue3","yellow","brown2","lightgreen","lightgoldenrod3","purple","pink","darkblue","gray","darkgreen"))( length(unique(tissues)) )[as.numeric(as.factor(tissues))]

		
		ann.col.new = matrix(ann.col,nr=length(thresholds),nc=length(ann.col),byrow=T)
		for (i in 1:length(thresholdsP)){		
			ann.col.new[i,which(DATA_p[i,]<tr)]=0
		}
		ord = order(tissues)
		col.thresh = colorRampPalette(c("black","firebrick3","tomato","RosyBrown2","dodgerblue3","skyblue2","lightskyblue1","gray70","blanchedalmond"))( length(unique(thresholds))+1 )
		if (length(thresholds)<4){
			rws = length(thresholds):1
		} else {
			rws = 4:1
		}
		tissues = gsub("_", " ", as.character(tissues))
		nms = gsub("_", " ", as.character(nms))
		pdf(paste(output_prefix,".",category,".pdf",sep=""),10,10)
		
		garfield.plot(DATA[,ord],ann.cols=ann.col.new[rws,ord], ann.col=ann.col[ord], ann.pch=15, rp.type="p",line.col=col.thresh,show.grid=TRUE, show.radial.grid=TRUE,labels=nms[ord],breaks=tissues[ord], radlab=TRUE,cex.axis=0.1, cex.lab=0.1, mar = c(2.5, 6.6, 3.5, 13.5), label.prop=1.1,poly.col=col.thresh, compact=compact)
		if (plot_title!=""){
			title(main=paste(plot_title," ",category,sep=""),line=7,cex=2) 
		}
		legend(-max(DATA,na.rm=T)*1.5,-max(DATA,na.rm=T)*1.6,c(thresholds,"1"),col=col.thresh,lty=1,lwd=6,title="P threshold",horiz=TRUE,cex=1,bty="n")
		dev.off()

	}

#### Output an exit message

print("Figures and table of results created!")
