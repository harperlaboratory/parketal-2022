# This workflow is used to quantify absolute amounts bythe  SIM method

# import =====
SkylineID <- read_csv('./import/SklineID_greek.csv') # peptide sequences and their ID





# PNS =====
PNS_fmol <- read_csv('./import/PNS_skyline.csv') # absolute amount of total 11 channels quantified from Skyline

PNS_ID_fmol <- inner_join(SkylineID, PNS_fmol, by = c('Skyline_ID' = 'Peptide'))

PNS_TMHQ <- read_csv('./import/PNS_isotope_norm_data_ID.csv') # relative quantitation from TOMAHAQ



PNS_ID_fmol_TMHQ <- inner_join(PNS_TMHQ, PNS_ID_fmol, by = 'ID')

PNS_ID_fmol_each <- PNS_ID_fmol_TMHQ[, 1] %>% 
  cbind(PNS_ID_fmol_TMHQ[, 2:13]/ PNS_ID_fmol_TMHQ$SumSN * PNS_ID_fmol_TMHQ$fmol) # absolute amount of each channel


PNS_ID_fmol_each_stat <- PNS_ID_fmol_each %>%
  cbind(DMSO_avg = PNS_ID_fmol_each[, 2:4] %>% rowMeans(),
        BSI_avg = PNS_ID_fmol_each[, 5:7] %>% rowMeans(),
        GSI_avg = PNS_ID_fmol_each[, 8:10] %>% rowMeans(),
        KO_avg = PNS_ID_fmol_each[, 11:12] %>% rowMeans()
        ) %>% 
  mutate('DMSO-KO' = DMSO_avg - KO_avg,
         'BSI-KO' = BSI_avg - KO_avg,
         'GSI-KO' = GSI_avg - KO_avg)





# PNS_Ami =====
PNS_Ami_fmol <- read_csv('./import/PNS_Ami_skyline.csv') # absolute amount of total 11 channels quantified from Skyline

PNS_Ami_ID_fmol <- inner_join(SkylineID, PNS_Ami_fmol, by = c('Skyline_ID' = 'Peptide'))



PNS_Ami_TMHQ <- read_csv('./import/PNS_Ami_isotope_norm_data_ID.csv') # relative quantitation from TOMAHAQ



PNS_Ami_ID_fmol_TMHQ <- inner_join(PNS_Ami_TMHQ, PNS_Ami_ID_fmol, by = 'ID')

PNS_Ami_ID_fmol_each <- PNS_Ami_ID_fmol_TMHQ[, 1] %>% 
  cbind(PNS_Ami_ID_fmol_TMHQ[, 2:13]/ PNS_Ami_ID_fmol_TMHQ$SumSN * PNS_Ami_ID_fmol_TMHQ$fmol) # absolute amount of each channel


PNS_Ami_ID_fmol_each_stat <- PNS_Ami_ID_fmol_each %>%
  cbind(DMSO_avg = PNS_Ami_ID_fmol_each[, 2:4] %>% rowMeans(),
        BSI_avg = PNS_Ami_ID_fmol_each[, 5:7] %>% rowMeans(),
        GSI_avg = PNS_Ami_ID_fmol_each[, 8:10] %>% rowMeans(),
        KO_avg = PNS_Ami_ID_fmol_each[, 11:12] %>% rowMeans()
        )  %>% 
  mutate('DMSO-KO' = DMSO_avg - KO_avg,
         'BSI-KO' = BSI_avg - KO_avg,
         'GSI-KO' = GSI_avg - KO_avg)




# Endo ===== 
Endo_fmol <- read_csv('./import/Endo_skyline.csv') # absolute amount of total 11 channels quantified from Skyline

Endo_ID_fmol <- inner_join(SkylineID, Endo_fmol, by = c('Skyline_ID' = 'Peptide')) 



Endo_TMHQ <- read_csv('./import/Endo_isotope_norm_data_ID.csv') # relative quantitation from TOMAHAQ


Endo_ID_fmol_TMHQ <- inner_join(Endo_TMHQ, Endo_ID_fmol, by = 'ID')

Endo_ID_fmol_each <- Endo_ID_fmol_TMHQ[, 1] %>% 
  cbind(Endo_ID_fmol_TMHQ[, 2:13]/ Endo_ID_fmol_TMHQ$SumSN * Endo_ID_fmol_TMHQ$fmol) # absolute amount of each channel


Endo_ID_fmol_each_stat <- Endo_ID_fmol_each %>%
  cbind(DMSO_avg = Endo_ID_fmol_each[, 2:4] %>% rowMeans(),
        BSI_avg = Endo_ID_fmol_each[, 5:7] %>% rowMeans(),
        GSI_avg = Endo_ID_fmol_each[, 8:10] %>% rowMeans(),
        KO_avg = Endo_ID_fmol_each[, 11:12] %>% rowMeans()
        ) %>% 
  mutate('DMSO-KO' = DMSO_avg - KO_avg,
         'BSI-KO' = BSI_avg - KO_avg,
         'GSI-KO' = GSI_avg - KO_avg)






# Endo_Ami =====
Endo_Ami_fmol <- read_csv('./import/Endo_Ami_skyline.csv') # absolute amount of total 11 channels quantified from Skyline

Endo_Ami_ID_fmol <- inner_join(SkylineID, Endo_Ami_fmol, by = c('Skyline_ID' = 'Peptide'))


Endo_Ami_TMHQ <- read_csv('./import/Endo_Ami_isotope_norm_data_ID.csv') # relative quantitation from TOMAHAQ



Endo_Ami_ID_fmol_TMHQ <-inner_join(Endo_Ami_TMHQ, Endo_Ami_ID_fmol, by = 'ID')
view(Endo_Ami_ID_fmol_TMHQ)

Endo_Ami_ID_fmol_each <- Endo_Ami_ID_fmol_TMHQ[, 1] %>% 
  cbind(Endo_Ami_ID_fmol_TMHQ[, 2:13]/ Endo_Ami_ID_fmol_TMHQ$SumSN * Endo_Ami_ID_fmol_TMHQ$fmol) # absolute amount of each channel


Endo_Ami_ID_fmol_each_stat <- Endo_Ami_ID_fmol_each %>%
  cbind(DMSO_avg = Endo_Ami_ID_fmol_each[, 2:4] %>% rowMeans(),
        BSI_avg = Endo_Ami_ID_fmol_each[, 5:7] %>% rowMeans(),
        GSI_avg = Endo_Ami_ID_fmol_each[, 8:10] %>% rowMeans(),
        KO_avg = Endo_Ami_ID_fmol_each[, 11:12] %>% rowMeans()
    ) %>% 
  mutate('DMSO-KO' = DMSO_avg - KO_avg,
         'BSI-KO' = BSI_avg - KO_avg,
         'GSI-KO' = GSI_avg - KO_avg)






# Lyso =====
Lyso_fmol <- read_csv('./import/Lyso_skyline.csv') # absolute amount of total 11 channels quantified from Skyline

Lyso_ID_fmol <- inner_join(SkylineID, Lyso_fmol, by = c('Skyline_ID' = 'Peptide'))



Lyso_TMHQ <- read_csv('./import/Lyso_isotope_norm_data_ID.csv') # relative quantitation from TOMAHAQ



Lyso_ID_fmol_TMHQ <- inner_join(Lyso_TMHQ, Lyso_ID_fmol, by = 'ID')

Lyso_ID_fmol_each <- Lyso_ID_fmol_TMHQ[, 1] %>% 
  cbind(Lyso_ID_fmol_TMHQ[, 2:13]/ Lyso_ID_fmol_TMHQ$SumSN * Lyso_ID_fmol_TMHQ$fmol) # absolute amount of each channel


Lyso_ID_fmol_each_stat <- Lyso_ID_fmol_each %>%
  cbind(DMSO_avg = Lyso_ID_fmol_each[, 2:4] %>% rowMeans(),
        BSI_avg = Lyso_ID_fmol_each[, 5:7] %>% rowMeans(),
        GSI_avg = Lyso_ID_fmol_each[, 8:10] %>% rowMeans(),
        KO_avg = Lyso_ID_fmol_each[, 11:12] %>% rowMeans()
    )  %>% 
  mutate('DMSO-KO' = DMSO_avg - KO_avg,
         'BSI-KO' = BSI_avg - KO_avg,
         'GSI-KO' = GSI_avg - KO_avg)






# Lyso_Ami =====
Lyso_Ami_fmol <- read_csv('./import/Lyso_Ami_skyline.csv') # absolute amount of total 11 channels quantified from Skyline

Lyso_Ami_ID_fmol <- inner_join(SkylineID, Lyso_Ami_fmol, by = c('Skyline_ID' = 'Peptide'))



Lyso_Ami_TMHQ <- read_csv('./import/Lyso_Ami_isotope_norm_data_ID.csv') # relative quantitation from TOMAHAQ



Lyso_Ami_ID_fmol_TMHQ <- inner_join(Lyso_Ami_TMHQ, Lyso_Ami_ID_fmol, by = 'ID')

Lyso_Ami_ID_fmol_each <- Lyso_Ami_ID_fmol_TMHQ[, 1] %>% 
  cbind(Lyso_Ami_ID_fmol_TMHQ[, 2:13]/ Lyso_Ami_ID_fmol_TMHQ$SumSN * Lyso_Ami_ID_fmol_TMHQ$fmol) # absolute amount of each channel


Lyso_Ami_ID_fmol_each_stat <- Lyso_Ami_ID_fmol_each %>%
  cbind(DMSO_avg = Lyso_Ami_ID_fmol_each[, 2:4] %>% rowMeans(),
        BSI_avg = Lyso_Ami_ID_fmol_each[, 5:7] %>% rowMeans(),
        GSI_avg = Lyso_Ami_ID_fmol_each[, 8:10] %>% rowMeans(),
        KO_avg = Lyso_Ami_ID_fmol_each[, 11:12] %>% rowMeans()
      )  %>% 
  mutate('DMSO-KO' = DMSO_avg - KO_avg,
         'BSI-KO' = BSI_avg - KO_avg,
         'GSI-KO' = GSI_avg - KO_avg)