library(DESeq2)
library(biomaRt)


# DE Analysis on Counts Data
path <- './GSE151879_raw_counts_genes.Macrophages.txt'
raw_data <- read.delim(path)

count_data <- raw_data[,2:7]
rownames(count_data) <- raw_data[,1] # gene names


# DDS data
condition <- c('Control', 'Control', 'Control', 'COVID_19', 'COVID_19', 'COVID_19') 
colData <- data.frame(row.names=colnames(count_data), condition=factor(condition, levels=c('Control','COVID_19')))

dds_data <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~condition)


# DESeq2 Analysis
dds <- DESeq(dds_data)
dds <- dds[!rowSums(counts(dds)) == 0, ] # exclude zero rows


# Results: adjusted p-value cutoff : 5%
result <- results(dds, contrast=c('condition','COVID_19','Control'), alpha = 0.05)

resultCleaned <- result[complete.cases(result$padj),]
# Significantly differentially expressed genes: Log2 fold changes and adjusted p-values are considered for this purpose. 
# Criterions:  |log2FoldChange| > 0  &  padj < 0.05  -->  4233 genes
sum(abs(resultCleaned$log2FoldChange)>1 & resultCleaned$padj < 0.05, na.rm = TRUE)


# Top 10 differentially expressed genes
ordered_result <- resultCleaned[order(resultCleaned$padj),]
top_genes <- ordered_result[abs(ordered_result$log2FoldChange) > 1 & ordered_result$padj < 0.05, ]
top10genes <- rownames(head(top_genes,10))

# Output:
# "ENSG00000004799"
# "ENSG00000150337"
# "ENSG00000162745"
# "ENSG00000149485"
# "ENSG00000184557"
# "ENSG00000074416"
# "ENSG00000134824"
# "ENSG00000111424"
# "ENSG00000103522"
# "ENSG00000059804"


# Finding gene names 
top_genes$ensembl <- rownames(top_genes)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl" )

# exon lengths in "bp", need "kbp" to normalize in TPM calc
# listAttributes(ensembl)

geneNames <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),  # only 1543 of 1607 found !
                  filters = "ensembl_gene_id",
                  values = top_genes$ensembl,
                  mart = ensembl )

idx <- match(top_genes$ensembl, geneNames$ensembl_gene_id)
top_genes$Names <- geneNames$hgnc_symbol[idx]
top_genes[c(top10genes), c('ensembl', 'Names')]

# Output:
# ENSG00000004799        PDK4
# ENSG00000150337      FCGR1A
# ENSG00000162745     OLFML2B
# ENSG00000149485       FADS1
# ENSG00000184557       SOCS3
# ENSG00000074416        MGLL
# ENSG00000134824       FADS2
# ENSG00000111424         VDR
# ENSG00000103522       IL21R
# ENSG00000059804      SLC2A3



# Finding Exon Lengths
exonLengths <- getBM( attributes = c("ensembl_gene_id", "exon_chrom_start", "exon_chrom_end"),
                      filters = "ensembl_gene_id",
                      values = rownames(count_data),
                      mart = ensembl )


# TPM calculation

# Multiple exon lengths on the genes
exonLengths$lengths <- (exonLengths[c('exon_chrom_end')] - exonLengths[c('exon_chrom_start')]) / 1000

# Total exon lengths wrt individual genes
tmp <- exonLengths[c('ensembl_gene_id', 'lengths')]
totalLengthsExon <- aggregate(tmp$lengths, by = list(tmp$ensembl_gene_id) , FUN = sum)

colnames(totalLengthsExon) <- c('GeneName', 'ExonLength') 


# RPK (reads per kilobase)
RPK <- count_data[c(totalLengthsExon$GeneName), ]/totalLengthsExon$ExonLength

# TPM
TPM_data <- RPK / colSums(RPK) * 1e+6

# TPM values to integer for DESeq2 analysis
for (col in colnames(TPM_data)) {
  TPM_data[[col]] = as.integer(TPM_data[[col]])
}



