# log2FC classified by lipid class ======
lipidomics_FC_class <- Endo_loadnorm_expnorm_expand %>% 
  ggplot(aes(y = log2FC, x = reorder(Class, -log2FC, median))) +  # x-axis: lipids are classified by 'Class', and the Classes are ordered by median of their elements
  geom_boxplot(width=0.7, lwd = 0.18, outlier.shape = NA) + # aesthetics of the boxplot
  theme_bw() +
  geom_jitter(width=0.2, aes(size = TSBH <= 0.01, color = TSBH <= 0.01, # 'TSBH' is the adjusted p-value. Adjust '0.01' to your target FDR
                             alpha = TSBH <= 0.01, shape = TSBH <= 0.01)) +
  scale_color_manual(values = c("TRUE" = "#BE6C29", "FALSE" = "gray20"), # change color codes here
                     name = expression(paste('Adjusted ', italic('p'), '-value')),
                     labels = c('> 0.01', expression(''<= 0.01)))+ 
  scale_alpha_manual(values = c(0.2, 0.8), name = NULL, labels = NULL, breaks = NULL) + # change transparency here
  scale_shape_manual(values = c(1, 16), name = NULL, labels = NULL, breaks = NULL) + # change shape here
  scale_size_manual(values = c(0.4, 0.6), name = NULL, labels = NULL, breaks = NULL) + # change size here
  xlab('Class') + # x-axis label
  scale_y_continuous(name =expression('Log'['2']*'(Endo-IP/Control-IP)'), + # y-axis label
                       limits = c(-5, 2.5), breaks = seq(-5, 2.5, by = 1), minor_breaks = NULL) + # y-axis scale
  geom_hline(yintercept= 0, col='black', size=0.15, alpha = 0.5) +
  theme(legend.position="bottom", text=element_text(size=8), axis.text = element_text(size =7),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(size = 0.2),
        axis.text.x = element_text(angle = 40, hjust=1))