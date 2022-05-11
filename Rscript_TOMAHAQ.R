# This workflow is used to analyze TOMAHAQ data

# peptide ID =======
peptideID <- read_csv('./import/peptideID_v2_greek.csv')

beta_peptide <- pull(peptideID[15:17, 2])
gamma_peptide <- pull(peptideID[20:26, 2])


# TMT isotopic impurity correction
isotope_correction <- read_csv('./import/TMT11plex_correction.csv')



# PNS =======
PNS_raw <- read_csv('./import/PNS_import.csv') #delete 1st col, insert 1st row with loading, annotate column name

PNS_isotope <- PNS_raw[, 1] %>% 
  cbind(as.matrix(PNS_raw[, 2:12]) %*% as.matrix(isotope_correction[, -1])) %>% 
  mutate(rowSums(.[, -1])) %>% 
  `colnames<-`(colnames(PNS_raw)) #isotope correction

PNS_normfactor <- 1 / as.numeric(t(PNS_isotope)[2:12, 1]) #inverse loading amount

PNS_isotope_norm <- PNS_isotope[, 1] %>% 
  cbind(as.data.frame(mapply('*', PNS_isotope[, 2:12], PNS_normfactor))) # inter-channel normalization

PNS_isotope_norm$. <- as.character(PNS_isotope_norm$.) #1st row to character

PNS_isotope_norm[PNS_isotope_norm < 0] <- 0 #replace minus to 0


PNS_isotope_norm_sumSN <- PNS_isotope_norm[, -1] %>% 
  rowSums() #sumSN

PNS_isotope_norm_data <- PNS_isotope_norm %>% #final table
  cbind(PNS_isotope_norm_sumSN) %>%
  `colnames<-`(colnames(PNS_raw)) %>%
  .[-1, ] #remove loading row



PNS_isotope_norm_data_ID <- 
  inner_join(peptideID, PNS_isotope_norm_data, by='Peptide') %>% select(-1) #peptide ID





# PNS_Ami =======
PNS_Ami_raw <- read_csv('./import/PNS_Ami_import.csv') #delete 1st col, insert 1st row with loading, annotate column name

PNS_Ami_isotope <- PNS_Ami_raw[, 1] %>% 
  cbind(as.matrix(PNS_Ami_raw[, 2:12]) %*% as.matrix(isotope_correction[, -1])) %>% 
  mutate(rowSums(.[, -1])) %>% 
  `colnames<-`(colnames(PNS_Ami_raw)) #isotope correction

PNS_Ami_normfactor <- 1 / as.numeric(t(PNS_Ami_isotope)[2:12, 1]) #inverse loading amount

PNS_Ami_isotope_norm <- PNS_Ami_isotope[, 1] %>% 
  cbind(as.data.frame(mapply('*', PNS_Ami_isotope[, 2:12], PNS_Ami_normfactor))) #channel norm

PNS_Ami_isotope_norm$. <- as.character(PNS_Ami_isotope_norm$.) #1st row to character

PNS_Ami_isotope_norm[PNS_Ami_isotope_norm < 0] <- 0 #replace minus to 0


PNS_Ami_isotope_norm_sumSN <- PNS_Ami_isotope_norm[, -1] %>% 
  rowSums() #sumSN

PNS_Ami_isotope_norm_data <- PNS_Ami_isotope_norm %>% #final table
  cbind(PNS_Ami_isotope_norm_sumSN) %>%
  `colnames<-`(colnames(PNS_Ami_raw)) %>%
  .[-1, ] #remove loading row


PNS_Ami_isotope_norm_data_ID <- 
  inner_join(peptideID, PNS_Ami_isotope_norm_data, by='Peptide') %>% select(-1) #peptide ID






# Endo =======
Endo_raw <- read_csv('./import/Endo_import.csv') #delete 1st col, insert 1st row with loading, annotate column name

Endo_isotope <- Endo_raw[, 1] %>% 
  cbind(as.matrix(Endo_raw[, 2:12]) %*% as.matrix(isotope_correction[, -1])) %>% 
  mutate(rowSums(.[, -1])) %>% 
  `colnames<-`(colnames(Endo_raw)) #isotope correction

Endo_normfactor <- 1 / as.numeric(t(Endo_isotope)[2:12, 1]) #inverse loading amount

Endo_isotope_norm <- Endo_isotope[, 1] %>% 
  cbind(as.data.frame(mapply('*', Endo_isotope[, 2:12], Endo_normfactor))) #channel norm

Endo_isotope_norm$. <- as.character(Endo_isotope_norm$.) #1st row to character

Endo_isotope_norm[Endo_isotope_norm < 0] <- 0 #replace minus to 0


Endo_isotope_norm_sumSN <- Endo_isotope_norm[, -1] %>% 
  rowSums() #sumSN

Endo_isotope_norm_data <- Endo_isotope_norm %>% #final table
  cbind(Endo_isotope_norm_sumSN) %>%
  `colnames<-`(c(colnames(Endo_raw), 'SumSN')) %>%
  .[-1, ] #remove loading row



Endo_isotope_norm_data_ID <- 
  inner_join(peptideID, Endo_isotope_norm_data, by='Peptide') %>% select(-1) #peptide ID





# Endo_Ami =======
Endo_Ami_raw <- read_csv('./import/Endo_Ami_import.csv') #delete 1st col, insert 1st row with loading, annotate column name

Endo_Ami_isotope <- Endo_Ami_raw[, 1] %>% 
  cbind(as.matrix(Endo_Ami_raw[, 2:12]) %*% as.matrix(isotope_correction[, -1])) %>% 
  mutate(rowSums(.[, -1])) %>% 
  `colnames<-`(colnames(Endo_Ami_raw)) #isotope correction

Endo_Ami_normfactor <- 1 / as.numeric(t(Endo_Ami_isotope)[2:12, 1]) #inverse loading amount

Endo_Ami_isotope_norm <- Endo_Ami_isotope[, 1] %>% 
  cbind(as.data.frame(mapply('*', Endo_Ami_isotope[, 2:12], Endo_Ami_normfactor))) #channel norm

Endo_Ami_isotope_norm$. <- as.character(Endo_Ami_isotope_norm$.) #1st row to character

Endo_Ami_isotope_norm[Endo_Ami_isotope_norm < 0] <- 0 #replace minus to 0


Endo_Ami_isotope_norm_sumSN <- Endo_Ami_isotope_norm[, -1] %>% 
  rowSums() #sumSN

Endo_Ami_isotope_norm_data <- Endo_Ami_isotope_norm %>% #final table
  cbind(Endo_Ami_isotope_norm_sumSN) %>%
  `colnames<-`(c(colnames(Endo_Ami_raw), 'SumSN')) %>%
  .[-1, ] #remove loading row



Endo_Ami_isotope_norm_data_ID <- 
  inner_join(peptideID, Endo_Ami_isotope_norm_data, by='Peptide') %>% select(-1) #peptide ID





# Lyso =======
Lyso_raw <- read_csv('./import/Lyso_import.csv') #delete 1st col, insert 1st row with loading, annotate column name

Lyso_isotope <- Lyso_raw[, 1] %>% 
  cbind(as.matrix(Lyso_raw[, 2:12]) %*% as.matrix(isotope_correction[, -1])) %>% 
  mutate(rowSums(.[, -1])) %>% 
  `colnames<-`(colnames(Lyso_raw)) #isotope correction

Lyso_normfactor <- 1 / as.numeric(t(Lyso_isotope)[2:12, 1]) #inverse loading amount

Lyso_isotope_norm <- Lyso_isotope[, 1] %>% 
  cbind(as.data.frame(mapply('*', Lyso_isotope[, 2:12], Lyso_normfactor))) #channel norm

Lyso_isotope_norm$. <- as.character(Lyso_isotope_norm$.) #1st row to character

Lyso_isotope_norm[Lyso_isotope_norm < 0] <- 0 #replace minus to 0


Lyso_isotope_norm_sumSN <- Lyso_isotope_norm[, -1] %>% 
  rowSums() #sumSN

Lyso_isotope_norm_data <- Lyso_isotope_norm %>% #final table
  cbind(Lyso_isotope_norm_sumSN) %>%
  `colnames<-`(colnames(Lyso_raw)) %>%
  .[-1, ] #remove loading row



Lyso_isotope_norm_data_ID <- 
  inner_join(peptideID, Lyso_isotope_norm_data, by='Peptide') %>% select(-1) #peptide ID






# Lyso_Ami =======
Lyso_Ami_raw <- read_csv('./import/Lyso_Ami_import.csv') #delete 1st col, insert 1st row with loading, annotate column name

Lyso_Ami_isotope <- Lyso_Ami_raw[, 1] %>% 
  cbind(as.matrix(Lyso_Ami_raw[, 2:12]) %*% as.matrix(isotope_correction[, -1])) %>% 
  mutate(rowSums(.[, -1])) %>% 
  `colnames<-`(colnames(Lyso_Ami_raw)) #isotope correction

Lyso_Ami_normfactor <- 1 / as.numeric(t(Lyso_Ami_isotope)[2:12, 1]) #inverse loading amount

Lyso_Ami_isotope_norm <- Lyso_Ami_isotope[, 1] %>% 
  cbind(as.data.frame(mapply('*', Lyso_Ami_isotope[, 2:12], Lyso_Ami_normfactor))) #channel norm

Lyso_Ami_isotope_norm$. <- as.character(Lyso_Ami_isotope_norm$.) #1st row to character

Lyso_Ami_isotope_norm[Lyso_Ami_isotope_norm < 0] <- 0 #replace minus to 0


Lyso_Ami_isotope_norm_sumSN <- Lyso_Ami_isotope_norm[, -1] %>% 
  rowSums() #sumSN

Lyso_Ami_isotope_norm_data <- Lyso_Ami_isotope_norm %>% #final table
  cbind(Lyso_Ami_isotope_norm_sumSN) %>%
  `colnames<-`(c(colnames(Lyso_Ami_raw), 'SumSN')) %>%
  .[-1, ] #remove loading row



Lyso_Ami_isotope_norm_data_ID <- 
  inner_join(peptideID, Lyso_Ami_isotope_norm_data, by='Peptide') %>% select(-1) #peptide ID