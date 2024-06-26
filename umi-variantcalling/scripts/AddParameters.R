#!/usr/bin/env Rscript
# Rscript AddParameters.R input.csv output.csv AML_candidates.txt AML_drivers.txt CHIP_mut.txt
args = commandArgs(trailingOnly=TRUE)

input <- read.table(args[1], 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
amlcand <- scan(file = args[3], 
                what = character())
amldriver <- scan(file = args[4], 
                  what = character())
chipmut <- read.table(args[5], 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
eb <- read.table(file = args[6], header = FALSE,
                 stringsAsFactors = FALSE)
eb <- as.vector(eb[,1])
eb[which(eb == ".")] <- NA



# Calculate Fisher exact test
# https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test
fisher <- round(dhyper(as.numeric(input["TR1_plus"][[1]]), 
                 as.numeric(input["TR1"][[1]]), 
                 as.numeric(input["TR2"][[1]]), 
                 as.numeric(input["TR1_plus"][[1]] + input["TR2_plus"][[1]]), 
                 log = TRUE), digits = 4) * -1
input["FisherScore"] <- fisher

input["MultiAllelic"] <- 0
input["EBScore"] <- 0

if (nrow(input) != length(eb))
{
    y <- which(colnames(input) == "MultiAllelic")
    z <- which(colnames(input) == "EBScore")
    input[1,z] <- eb[1]
    j = 2
    for (i in 2:nrow(input))
    {
        if (all(input[i-1,c(1:3)] == input[i,c(1:3)]))
        {
            input[i,z] <- eb[j]
            input[i-1,y] <- 1
            input[i,y] <- 1
            
        } else {
            input[i,z] <- eb[j]
            j = j + 1
        }
    }
    
    
} else {
    input["EBScore"] <- eb
    
}



# Calculate VAF
input["TVAF"] <- round(as.numeric(input["TR2"][[1]])/
                           as.numeric(input["readDepth"][[1]]), digits = 4)

icand <- vector(length = nrow(input))
idriver <- vector(length = nrow(input))

for (i in 1:nrow(input))
{
    icand[i] <- any(strsplit(input["Gene.refGene"][[1]][i], split = ";")[[1]] %in% amlcand)
    idriver[i] <- any(strsplit(input["Gene.refGene"][[1]][i], split = ";")[[1]] %in% amldriver)
}

input["AMLcandidate"] <- icand
input["AMLdriver"] <- idriver

colnames(input)[which(colnames(input) == "Gene.refGene")] <- "Gene"
colnames(input)[which(colnames(input) == "GeneDetail.refGene")] <- "GeneDetail"
colnames(input)[which(colnames(input) == "ExonicFunc.refGene")] <- "ExonicFunc"
colnames(input)[which(colnames(input) == "AAChange.refGene")] <- "AAChange"

output <- merge(input, chipmut, all.x = TRUE)[union(names(input), names(chipmut))]

write.table(output, file = args[2], quote = FALSE, 
            row.names = FALSE, sep = "\t")