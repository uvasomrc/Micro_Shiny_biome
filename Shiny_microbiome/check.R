countdf <- read.table("count_data_relman2017_samples.otu_table.txt", header=TRUE,sep = "\t" )
taxodf <- read.table("taxonomy_relman2017_samples.tax_table.txt", header=TRUE,sep = "\t" )
sampledf <- read.table("sampe_data_relman2017_samples.sample_data.txt", header=TRUE,sep = "\t" )

nrow(countdf) == nrow(taxodf)

ncol(countdf) == nrow(sampledf)
