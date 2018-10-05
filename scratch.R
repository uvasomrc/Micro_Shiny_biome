library(readr)
library(dplyr)
library(tibble)
library(reshape2)


#library(stringr)

cdf <- read.table("relman2017_samples.otu_table.txt", header=TRUE,sep = "\t" )
tdf <- read.table("relman2017_samples.tax_table.txt", header=TRUE,sep = "\t" )
sdf <- read.table("relman2017_samples.sample_data.txt", header=TRUE,sep = "\t" )
cdf2 <- cdf[order(cdf$OTUId),]
rownames(cdf2) <- cdf2[,1]
cdf2[,1] <- NULL
tdf2 <- tdf[order(tdf$X),]
rownames(tdf2) <- tdf2[,1]
tdf2[,1] <- NULL

sdf2 <- sdf[order(sdf$sample),]
rownames(sdf2) <- sdf2[,1]
sdf2[,1] <- NULL



table(tdf2$domain)
table(tdf2$phylum)
table(tdf2$class)
table(tdf2$order)
table(tdf2$family)
table(tdf2$genus)


#######

# taxo_raw <- as.data.frame(original_data()$tdf_orig[input$taxo_level])
# taxo_raw <- as.data.frame(lapply(taxo_raw , function(x) {gsub(".*unclassified.*", "Unclassified", x)}))
# cdf_orig_raw <- original_data()$cdf_orig
# cdf_orig_raw <- as.data.frame(lapply(cdf_orig_raw, as.numeric))
# df_raw <- cbind(cdf_orig_raw, taxo_raw)
# 
# df_raw2 <- df_raw %>% 
#   
#   group_by(input$taxo_level) %>%
#   summarise_all(funs(sum)) %>%
#   as.data.frame()
# 
# df_raw2


xtdf2 <- as.data.frame(tdf2$class)
colnames(xtdf2) <- c("taxo_class")
xtdf2 <-  as.data.frame(lapply(xtdf2, function(x) {gsub(".*unclassified.*", "Unclassified", x)}))
xcdf2 <- cdf2
xcdf2 <- as.data.frame(lapply(xcdf2, as.numeric))
xdf1 <- cbind(xcdf2, xtdf2)
xdf2 <- xdf1 %>%
  group_by(taxo_class) %>% 
  summarise_all(funs(sum)) %>% 
  as.data.frame()
rownames(xdf2) <- xdf2[,1]
xdf2[,1] <- NULL


df2_m <- melt(cbind(df2, ind=rownames(df2)), id.vars = c('ind'))

ggplot(df2_m,aes(x = variable, y = value, fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") +
  scale_y_continuous(labels = percent_format())
###########

percent_format