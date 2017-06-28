#script for handling Timema transcriptome polymorphisms

#get data

#SNPs
SNPs <- read.table (file="snp_per_gene_counts.csv", sep=",", header=TRUE)
#degenerate positions
degpos <- read.table (file="Tall_degenerate_pos.csv", sep="\t", header=TRUE)

#only select the degpos that are in SNPs
degpos <- degpos[degpos$transcript_id %in% SNPs$transcript, ] 

#convert to numeric again
degpos$f1 <- as.numeric(levels(degpos$f1))[degpos$f1]
degpos$f2 <- as.numeric(levels(degpos$f2))[degpos$f2]
degpos$f3 <- as.numeric(levels(degpos$f3))[degpos$f3]
degpos$f4 <- as.numeric(levels(degpos$f4))[degpos$f4]

#calculate syn and nonsyn sites
degpos$syn <- (degpos$f4+((degpos$f2+degpos$f3)/2))
degpos$nonsyn <- (degpos$f1+((degpos$f2+degpos$f3)/2))

#merge with snp data
polym <- as.data.frame(c(SNPs, degpos[c(2:7)]))

#calculate pN pS and pN/pS
polym$pS <- (polym$syn_total/polym$syn)
polym$pN <- (polym$ns_total/polym$nonsyn)

polym$pNpS <- (polym$pN/polym$pS)
