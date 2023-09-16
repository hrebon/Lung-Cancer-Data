library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)

# get a list of projects

gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LUAD')

## building a query
queryTCGA <- GDCquery(project = 'TCGA-LUAD',
                      data.category = 'Transcriptome Profiling')

#output_query_TCGA <- getResults(queryTCGA)

# build a query a retrieve gene expression data

queryTCGA <- GDCquery(project = 'TCGA-LUAD',
                      data.category = 'Transcriptome Profiling',
                      data.type = 'Gene Expression Quantification',
                      experimental.strategy = 'RNA-Seq',
                      workflow.type = 'STAR - Counts',
                      access = 'open')

getResults(queryTCGA)
# downlaod gene expression data of LUAD--
GDCdownload(queryTCGA)
# prepare data
tcga_luad_data <- GDCprepare(queryTCGA, summarizedExperiment = TRUE)
output <- assay(tcga_luad_data)

#downlaod DNA methylation data for lung cancer(LUAD)---

query_methly <- GDCquery(project = 'TCGA-LUAD',
                         data.category = 'DNA Methylation',
                         data.type = 'Methylation Beta Value',
                         platform = 'Illumina Human Methylation 450',
                         access = 'open')

output_query_methyl <- getResults(query_methly)


GDCdownload(query_methly)

# download mutation data of lung cancer

query_mutation <- GDCquery(project = 'TCGA-LUAD',
                           data.category = 'Simple Nucleotide Variation',
                           access = 'open',
                           )
getResults(query_mutation)

GDCdownload(query_mutation)

# download clinical data of lung cancer (LUAD)

query_clinical <- GDCquery(project = 'TCGA-LUAD',
                           data.category = 'Clinical',
                           data.type = 'Clinical Supplement',
                           data.format = 'xml',
                           access = 'open')

getResults(query_clinical)

GDCdownload(query_clinical)
