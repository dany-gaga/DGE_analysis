#Load packages----
library(limma)
library(edgeR)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#set the working directory
setwd('/home/dany/Desktop/G-Avenir/data-analysis')
# load the data
x<-read.delim("featurecount.csv",skip=0, sep=",", check.names=FALSE, row.names = 1)
# you can add the name of your own column names
counts <- x[,c('BD01_S9','BD02_S8','BD03_S7','BP02_S9','BP03_S8','BP04_S7','BU01_S6','BU02_S5','BU03_S4','DA01_S9','DA02_S8','DA03_S7','DD01_S6','DD02_S5','DD03_S4','DP01_S6','DP02_S5','DP03_S4','DU01_S3','DU02_S2','DU03_S1','KI01_S3','KI02_S2','KI03_S1','KS6_S3','KS9_S2','KS10_S1')]

# create list objects
list.a <- as.list(c("BD","BP","BU","DA","DD","DP","DU","KI","KS"))
list.b <- as.list(c("BD","BP","BU","DA","DD","DP","DU","KI","KS"))

# create an empty list to populate     
list.ab <- list()

# double loop to create recursive combinations using append 
    ct=0    
      for( i in 1:length(list.a) ) {    
        for (j in 1:length(list.b) ) {
          ct=ct+1
           list.ab[[ct]] <- append(list.a[[i]], list.b[[j]])     
         } 
       }
#create a new empty list
list.new <- list()
c=0
# remove all duplicate pais combinations from list.ab
for (col in list.ab) {
	if (col[1] != col[2]) {
		c=c+1
		list.new[[c]] <- append(col[1], col[2])
	}
}

# for all pairs combinations in the list.new, create a summary table, a maplot and a volcano plot
for (cn in list.new) {
	#create a subset table containing columns of the name in col and save them in a csv file
	sub_tab <- counts %>% dplyr:: select(grep(cn[1], names(counts)), grep(cn[2], names(counts))) 
	# filtering out all genes with max count < 50
	keep <- apply(counts, 1, max) >= 50
	x <- x[keep,]
	sub_tab <- sub_tab[keep,]
	csv_name <- paste0(cn[1], sep="_", cn[2], ".csv")
	write.csv(sub_tab, csv_name)
	# Defining a design matrix based on the experimental design.
	design <- matrix(data=c(1,1,1,0,0,0,0,0,0,1,1,1), nrow=6, ncol=2, dimnames = list(c(), c(cn[1],cn[2])))
	#creating a matrix of contrasts, where which each column represents a contrast between two groups of interest.
	cont.matrix <- matrix(data=c(-1,1), nrow=2, ncol=1, dimnames = list(c(), c(cn[1])))
	# Create a DGEList object
	y <- DGEList(sub_tab)
	#factor normalization
	y <- calcNormFactors(y, method="TMM")
	#estimate dispersion
	y <- estimateDisp(y,design)
	# Now conduct glm fit test for treatment effect
	fit <- glmFit(y,design)
	lrt <- glmLRT(fit, contrast=cont.matrix)
	#The total number of genes significantly up-regulated or down-regulated at 5% FDR is summarized as follows:
	# The number of down and up-regulated genes based on different FDR and logFC  could be also computed  
	c <- summary(decideTests(lrt, p.value = 0.05, lfc=0))
	csv_name1 <- paste0("summary", sep="_", cn[1], "vs", cn[2], ".csv")
	write.csv(c, csv_name1)
	c1 <- summary(decideTests(lrt, p.value = 0.01, lfc=0))
	csv_name2 <- paste0("summary1", sep="_", cn[1], "vs", cn[2], ".csv")
	write.csv(c1, csv_name2)
	c2 <- summary(decideTests(lrt, p.value = 0.05, lfc=2))
	csv_name3 <- paste0("summary2", sep="_", cn[1], "vs", cn[2], ".csv")
	write.csv(c2, csv_name3)
	c3 <- summary(decideTests(lrt, p.value = 0.01, lfc=2))
	csv_name4 <- paste0("summary3", sep="_", cn[1], "vs", cn[2], ".csv")
	write.csv(c3, csv_name4)
	# the top ten  differential expressed genes
	topTags(lrt)
	FDR <- p.adjust(lrt$table$PValue, method="BH")
	sum(FDR < 0.05)
	#Generate the name to be use to save the ma_plot as a png
	plot_name1 <- paste0("ma", sep="_", cn[1], "vs", cn[2], ".png")
	plot_name2 <- paste0(cn[1], sep="_", cn[2])
	#Generate the MA plot and save it as a png file
	png(file=plot_name1)
	plotMD(lrt, main=plot_name2)
	abline(h=c(-1,1), col="blue")
	dev.off()
	#generate the result table and save it save it in a file
	out <- topTags(lrt, n=Inf, sort.by='PValue')$table
	#save the name of the result table to be use
	csv_name2 <- paste0("stat1", sep="_", cn[1], sep="_", cn[2], ".csv")
	write.csv(out, csv_name2)
	#Generate the name to be use to save the volcano_plot as a png
	plot_name3 <- paste0("volcano", sep="_", cn[1], "vs", cn[2], ".png")
	#Generate the volcano plot
	d1 <- out
	# add a column of NAs
	d1$Legend <- "NotSig"
	# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
	d1$Legend[d1$logFC >= 1 & d1$PValue < 0.05] <- "Up"
	# if log2Foldchange < -1 and pvalue < 0.05, set as "Down"
	d1$Legend[d1$logFC <= -1 & d1$PValue < 0.05] <- "Down"
	# name of the Volcano plot
	plot_name4 <- paste0(cn[1], sep="_", cn[2])
	# Re-plot but this time color the points with "diff. expressed"
	ggplot(data=d1, aes(x=logFC, y=-log10(PValue), col=Legend)) + geom_point() +
	  theme(panel.background = element_rect(fill = "white", color = "grey50")) +
	  scale_color_manual(breaks = c("Down", "NotSig", "Up"),values=c("red", "black", "blue")) + 
	  geom_vline(xintercept=c(-1, 1), col="red", lty=3) + geom_hline(yintercept = -log10(0.05), lty=3, color="red") +
	  ggtitle(plot_name4)
	ggsave(file=plot_name3)
}	
	


