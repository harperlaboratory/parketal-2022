# This workflow is used to analyze whole cell, Endo-IP, and Lyso-IP proteomics

# import data =========
raw_table <- read_tsv('./import/protein_quant_25450.tsv')


# data clean-up ====

# reverse hits (##)
grep('^##', raw_table$`Protein Id`, value = T)
reverse_hit <- grep('^##', raw_table$`Protein Id`)


# contaminants (_contaminant)
grep('contaminant', raw_table$`Protein Id`, value = T)
contaminant <- grep('contaminant', raw_table$`Protein Id`)

data_filtered <- raw_table[-c(reverse_hit, contaminant), ]


# filter out NA
sum(is.na((data_filtered[,1])))


# extract protein identifier (._HUMAN) & sum
uniprotID <- strsplit(data_filtered$`Protein Id`, '\\|')
uniprotID_split <- matrix(unlist(uniprotID), ncol=3, byrow=T)
protID <- matrix(unlist(strsplit(uniprotID_split[,3], '_')), ncol=2, byrow=T)
ID <- protID[,1]



# Uniprot ID to gene name
ID_quant_sum <- uniprotID_split %>% 
  cbind(data_filtered[, c(2,3, 5, 8:15)]) %>% 
  `colnames<-`(c('Database', 'Entry', 'Entry name',
                 'Gene', 'Description', 'Number of peptides',
                 'FLAG_Ctrl_1', 'FLAG_Ctrl_2', 'FLAG_Ctrl_3', 'FLAG_Ctrl_4',
                 'FLAG_EEA1_1', 'FLAG_EEA1_2', 'FLAG_EEA1_3', 'FLAG_EEA1_4')) %>%
  .[, -1]




## low SN
summary(rowSums(ID_quant_sum[, 6:13]))

SNcutoff <- 100
SNsum <- rowSums(ID_quant_sum[, 6:13])
sum(SNsum < SNcutoff)


ID_quant_sum_SNfilter <- ID_quant_sum %>% mutate(SNsum) %>% 
  filter(SNsum >= SNcutoff) %>% select(-SNsum)




# channel normalization ====
loading_level <- colSums(ID_quant_sum_SNfilter[, 6:13])

loading_ratio <- loading_level / max(loading_level)

norm_factor <- 1/loading_ratio

as.vector(norm_factor)



normdata <- as.data.frame(ID_quant_sum_SNfilter[, 1:5])


for (i in 1:8) {
  normdata[, i+5] <- as.data.frame(ID_quant_sum_SNfilter[, i+5] * as.vector(norm_factor)[i])
}

names(normdata) <- names(ID_quant_sum)



# log2 transformation ====
normdata_log2 <- ID_quant_sum_SNfilter[, 1:5] %>% 
  cbind(log2(ID_quant_sum_SNfilter[, 6:13]))


## remove -Inf
sum(is.infinite(rowSums(normdata_log2[, -c(1:5)])))
normdata_log2_noInf <- normdata_log2[!is.infinite(rowSums(normdata_log2[, -c(1:5)])),]




# Itzhak proteome ===============
ItzhakProteome <- read_csv('./import/Itzhak_proteome.csv')
View(ItzhakProteome)


ItzhakProteome_confident <- ItzhakProteome %>% 
  filter(Confidence == 'Very High' | Confidence == 'High') %>% 
  select(Entry, Compartment)


normdata_log2_noInf_compartment <- normdata_log2_noInf %>% 
  left_join(ItzhakProteome_confident, by = c('Entry' = 'Entry')) %>% 
  replace_na(list(Compartment = "Nonspecific"))




# Endo IP analysis ====
  ## long format

Endo_data <- normdata_log2_noInf_compartment

Endo_data_long <- melt(data = Endo_data[, c(1, 6:13)],
                       id.vars = 'Entry', 
                       variable.name = "sample", 
                       value.name = "log2")

View(Endo_data_long)

Endo_data_group <- data.frame('sample'=c('FLAG_Ctrl_1', 'FLAG_Ctrl_2', 'FLAG_Ctrl_3', 'FLAG_Ctrl_4',
                                         'FLAG_EEA1_1', 'FLAG_EEA1_2', 'FLAG_EEA1_3', 'FLAG_EEA1_4'), 
                              'group'=as.factor(c(rep('FLAG_Ctrl', 4), rep('FLAG_EEA1', 4))))

Endo_data_long_group <- left_join(Endo_data_long, Endo_data_group, by='sample')



  ## multiple t-test

Endo_data_test <- Endo_data_long_group %>% 
  group_by(Entry) %>% 
  t_test(formula = log2 ~ group, var.equal =T) # Student's T-test



  # multtest package, HKY (2006) FDR method for Endo_data

Endo_data_test_TSBHadjp <- mt.rawp2adjp(Endo_data_test$p, proc=("TSBH"), 
                                        alpha = 0.05, na.rm = FALSE)

Endo_data_test_TSBH <- data.frame('Entry'=Endo_data_test$Entry, 
                                  Endo_data_test_TSBHadjp$adj[order(Endo_data_test_TSBHadjp$index), ])

Endo_data$Entry <- as.character(Endo_data$Entry)


Endo_data_expand <- Endo_data %>% 
  mutate(FLAG_Ctrl = rowMeans(Endo_data[, 6:9])) %>% 
  mutate(FLAG_EEA1 = rowMeans(Endo_data[, 10:13])) %>% 
  mutate(log2FC = FLAG_EEA1 - FLAG_Ctrl) 


Endo_data_test_TSBH$Entry <- as.character(Endo_data_test_TSBH$Entry)

Endo_data_expand_TSBH <- Endo_data_expand %>% 
  left_join(Endo_data_test_TSBH[, c(1,3)], by='Entry')


Endo_data_expand_TSBH <- Endo_data_expand_TSBH %>% replace_na(list(Compartment = "Nonspecific"))
