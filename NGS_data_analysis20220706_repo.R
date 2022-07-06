######################################################################

wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("C:/Users/ninag/Documents/Project_UTI-in-RTP")
library(tidyverse)
library(ggpubr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("lefser")

library(lefser)

### saving images
today.date <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d")
dir.create(paste0("R.export/plots/", today.date), showWarnings = FALSE, recursive = TRUE)
image.save.path <- paste0(wd, "/", paste0("R.export/plots/", today.date))

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), strip.background = element_rect(fill = "hotpink"), legend.background = element_blank(), legend.key = element_blank()) #text=element_text(family = "Courier")

openxlsx::getSheetNames("/Users/ninag/Documents/Project_UTI-in-RTP")
openxlsx::getSheetNames("/Users/ninag/Documents/Project_UTI-in-RTP/UTI_16s_NGS_data.xlsx")
data <- openxlsx::read.xlsx("/Users/ninag/Documents/Project_UTI-in-RTP/UTI_16s_NGS_data.xlsx", sheet = "geni_UTI")

str(data)

data2 <-
  data %>%
  dplyr::filter(rank == "genus")



# do.call(rbind, ...) instead of dplyr::bind_rows, which does not catch different column names though
lineage <- dplyr::bind_rows(lapply(strsplit(gsub("\"", "", data2$lineage), ";"), function(x) {
  temp <- rev(x)

  genus <- ifelse("genus" %in% temp, temp[which(temp == "genus")+1], NA)
  family <- ifelse("family" %in% temp, temp[which(temp == "family")+1], NA)
  order = ifelse("order" %in% temp, temp[which(temp == "order")+1], NA)
  class = ifelse("class" %in% temp, temp[which(temp == "class")+1], NA)
  phylum = ifelse("phylum" %in% temp, temp[which(temp == "phylum")+1], NA)

  df <- data.frame("genus" = genus, "family" = family,
                   "order" = order, "class" = class,
                   "phylum" = phylum, stringsAsFactors = F)
  return(df)
}))



### how to calculate relative abundance and how to pay respect to different sequencing depths is not trivial, or at least controversial
data3 <-
  data2 %>%
  dplyr::select(-lineage) %>%
  dplyr::bind_cols(lineage) %>%
  tidyr::pivot_longer(cols = dplyr::contains("UTI"), names_to = "sample", values_to = ("abundance")) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(total_abundance_per_sample = sum(abundance)) %>%
  dplyr::mutate(rel_genus_abundance = abundance/total_abundance_per_sample)
## aggregated values (relative abundance) in data3 refers to genus!

## do it by genus:
data4 <-
  data3 %>%
  dplyr::group_by(sample, genus, total_abundance_per_sample) %>%
  dplyr::summarise(abundance_per_genus = sum(abundance), .groups = "drop") %>%
  dplyr::mutate(rel_genus_abundance = abundance_per_genus/total_abundance_per_sample)
# order phylum (make factor) to have beautiful colors for most abundant phyla
#data4$phylum <- factor(data4$phylum, levels = data4 %>% dplyr::group_by(phylum) %>% dplyr::summarise(overall_abundance = sum(abundance_per_phylum)) %>% dplyr::arrange(-overall_abundance) %>% dplyr::pull(phylum))
str(data4)

#############################################
## ---- meta_data_join --------
sd <- openxlsx::read.xlsx("UTI_tidy-table_20220329.xlsx", detectDates = TRUE)
sd <- sd[,colSums(is.na(sd))<nrow(sd)]
ps <- dplyr::left_join(data4, sd, by = c("sample" = "ID"))

##########USE Lefser

table(ps$sex)

table(ps$sex, ps$Event)

res <- lefser(ps, groupCol = "Event", blockCol = "sex")

head(res)

lefserPlot(res)