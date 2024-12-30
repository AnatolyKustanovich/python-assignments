setwd("/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA")
getwd()

library(dplyr)
library(ggplot2)

#Read files
BACSinput = readRDS("794_MW293TfullMaximaShishkinBACSInput_S12_Aligned.sortedByCoord.out.txDT.rds")
BACStreat=readRDS("794_MW293TfullMaximaShishkinBACSTreat_S11_Aligned.sortedByCoord.out.txDT.rds")
#Change column "-" to "Del"
colnames(BACSinput)[colnames(BACSinput) == "-"] <- "Del"
colnames(BACStreat)[colnames(BACStreat) == "-"] <- "Del"

#Frequencies of nucleotides for treat file
A_freq <- c()
C_freq <- c()
G_freq <- c()
T_freq <- c()
C_conv <- c()
bkg_C_conv <- c()
bkg_G_conv <- c()
bkg_T_conv <- c()
bkg_A_conv <- c()
for (row in 1:nrow(BACStreat)) {
  cov1 = BACStreat$A[row] + BACStreat$G[row] + BACStreat$C[row] + BACStreat$T[row] 
  A_ratio <- round(BACStreat$A[row] / cov1 *100, digits = 1)
  A_freq <- append(A_freq, A_ratio)
  C_ratio <- round(BACStreat$C[row] / cov1*100, digits = 1)
  C_freq <- append(C_freq, C_ratio)
  G_ratio <- round(BACStreat$G[row] / cov1*100, digits = 1)
  G_freq <- append(G_freq, G_ratio)
  T_ratio <- round(BACStreat$T[row] / cov1*100, digits = 1)
  T_freq <- append(T_freq, T_ratio)
}
for (row in 1:nrow(BACStreat)) {
  cov3 = BACStreat$C[row] + BACStreat$T[row] 
  A_ratio <- round(BACStreat$A[row] / cov3 *100, digits = 1)
  bkg_A_conv <- append(bkg_A_conv, A_ratio)
  C_ratio <- round(BACStreat$C[row] / cov3 *100, digits = 1)
  bkg_C_conv <- append(bkg_C_conv, C_ratio)
  G_ratio <- round(BACStreat$G[row] / cov3 *100, digits = 1)
  bkg_G_conv <- append(bkg_G_conv, G_ratio)
  T_ratio <- round(BACStreat$T[row] / cov3 *100, digits = 1)
  bkg_T_conv <- append(bkg_T_conv, T_ratio)
}
for (row in 1:nrow(BACStreat)) {
  cov2 = BACStreat$C[row] + BACStreat$T[row] 
  C_convrat <- round(BACStreat$C[row] / cov2*100, digits = 1)
  C_conv <- append(C_conv, C_convrat)
}
#Add columns with frequencies
BACStreat$A_freq = A_freq
BACStreat$C_freq = C_freq
BACStreat$G_freq = G_freq
BACStreat$T_freq = T_freq
BACStreat$C_conv = C_conv
BACStreat$bkg_A_conv = bkg_A_conv
BACStreat$bkg_C_conv = bkg_C_conv
BACStreat$bkg_G_conv = bkg_G_conv
BACStreat$bkg_T_conv = bkg_T_conv

#Work with Taoka dataset
taoka = as.data.frame(read.table('human_Taoka_Modifications1.bed',header = FALSE, sep="\t",stringsAsFactors=FALSE))

colnames(taoka)[colnames(taoka) == "V2"] <- "gencoor"
colnames(taoka)[colnames(taoka) == "V5"] <- "refSeq"
colnames(taoka)[colnames(taoka) == "V1"] <- "chr"
colnames(taoka)[colnames(taoka) == "V4"] <- "Conv_Taoka"

BACStreat$common <- paste0(BACStreat$gene, "_", BACStreat$gencoor)
taoka$common <- paste0(taoka$chr, "_", taoka$gencoor)
taoka$common <- gsub("S", "s", taoka$common)
taoka$chr <- gsub("S", "s", taoka$chr)

#Merge 
BACStreat_taoka = merge(x=taoka, y=BACStreat, by="common", all.x = T, all.y = T )
write.csv(BACStreat_taoka,"~/BACStreat_taoka.csv", row.names = FALSE)

ggplot(BACStreat_taoka %>% filter(V3 == "Y", gene == "5.8s"), aes(x = gencoor.y)) +
  # Add segments for C_conv and Conv_Taoka
  geom_segment(aes(xend = gencoor.y, y = 0, yend = Conv_Taoka), color = "darkviolet") +
  geom_segment(aes(xend = gencoor.y, y = 0, yend = C_conv), color = "brown2") +
  
  # Add points for Conv_Taoka
  geom_point(aes(y = Conv_Taoka), color = "darkviolet", size = 3) +
  
  # Add points for C_conv
  geom_point(aes(y = C_conv), color = "brown2", size = 3) +
  
    # Add text labels for gencoor.x (optional for clarity)
  #geom_text(aes(y = C_conv, label = paste0(gencoor.x)), 
   #         vjust = -0.5, color = "black", size = 2) +
  #geom_text(aes(y = Conv_Taoka, label = paste0(gencoor.x)), 
   #         vjust = -0.5, color = "black", size = 2) +
  # Adjust theme and labels
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position in 5.8s") +
  ylab("Conversion rate (%)") +
  ggtitle("BACS. Comparison of C_conv and Conv_Taoka Rates")


ggplot(BACStreat_taoka %>% filter(refSeq.y == "T", gene == "28s"), aes(x = gencoor.y, shape = V3 == "Y")) +
  # Add segments for C_conv and Conv_Taoka
  geom_segment(aes(xend = gencoor.y, y = 0, yend = Conv_Taoka), color = "darksalmon") +
  geom_segment(aes(xend = gencoor.y, y = 0, yend = C_conv), color = "brown2") +
  
  geom_point(aes(y = Conv_Taoka), color = "darksalmon", size = 3) +
  
  geom_point(aes(y = C_conv), color = "brown2", size = 3) +
  
  # Define custom shapes for logical values
  scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 16)) + # Square for FALSE, Circle for TRUE
  
  # Adjust theme and labels
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position in 28s") +
  ylab("Conversion rate (%)") +
  ggtitle("BACS: Comparison of C_conv and Conv_Taoka Rates")


ggplot(BACStreat_taoka %>% filter(refSeq.y == "T", gene == "18s"), aes(x = gencoor.y)) +
  # Add segments for C_conv and Conv_Taoka
  geom_segment(aes(xend = gencoor.y, y = 0, yend = Conv_Taoka), color = "darkviolet") +
  geom_segment(aes(xend = gencoor.y, y = 0, yend = C_conv), color = "brown2") +
  
  # Add points for Conv_Taoka
  geom_point(aes(y = Conv_Taoka, shape = V3 == "Y"), color = "darkviolet", size = 2) +
  
  # Add points for C_conv
  geom_point(aes(y = C_conv, shape = V3 == "Y"), color = "brown2", size = 2) +
  
  # Add text labels for gencoor.x (optional for clarity)
  geom_text(aes(y = C_conv, label = paste0(gencoor.x)), 
           vjust = -0.5, color = "black", size = 2) +
  geom_text(aes(y = Conv_Taoka, label = paste0(gencoor.x)), 
           vjust = -0.5, color = "black", size = 2) +
  scale_shape_manual(values = c(FALSE = 15, TRUE = 16)) + # Square for FALSE, Circle for TRUE
  # Adjust theme and labels
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position in 18s") +
  ylab("Conversion rate (%)") +
  ggtitle("BACS. Comparison of C_conv and Conv_Taoka Rates")



colnames(treat_taoka)

library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)



# Plot the data
ggplot(treat_T_taoka%>%filter(V3== "Y"), aes(x = gencoor.y, y = bkg_C_conv, color = V3)) +
  geom_line() +
  theme_light() +
  xlab("Genomic Coordinates") +
  ylab("Background Conversion Rate") +
  ggtitle("Background Conversion Rates for T Reference Sequence")

#Background plot 
ggplot(treat_taoka%>% filter(refSeq.y == "T", chr.y == "18s"), aes(x = gencoor.y, y = bkg_C_conv, color = V3)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("Genomic Coordinates") +
  ylab("Conversion Rate (%)") +
  ggtitle("BACS. Conversion Rates for T Reference Sequence. 18s")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  ) +
  scale_color_manual(
    values = c(
      "m1acp3Y" = "purple", 
      "pU" = "blue", 
      "Um" = "green", 
      "Y" = "red", 
      "NA" = "black"
    )
  )+
  ylim(0, 100)

ggplot(treat_18s_taoka%>% filter(refSeq.y == "T"), aes(x = gencoor, y = bkg_A_conv, color = V3)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("Genomic Coordinates") +
  ylab("Conversion Rate (%)") +
  ggtitle("Conversion Rates for T Reference Sequence")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  ) +
  scale_color_manual(
    values = c(
      "m1acp3Y" = "red", 
      "pU" = "blue", 
      "Um" = "green", 
      "Y" = "purple", 
      "NA" = "black"
    )
  )+
ylim(0, 100)

  
ggplot(treat_18s%>% filter(refSeq == "T"), aes(x =gencoor, y = bkg_T_conv, color = variable)) +
  geom_line()

plot <- ggplot(treat_18s_T, aes(x=gencoor)) +
  #Three geom_segments 
  geom_segment(aes(y=C_conv,yend=C_conv,xend=gencoor), colour= "blue") +
  geom_segment(aes(y=C_freq,yend=C_freq,xend=gencoor), colour= "red") +
  geom_segment(aes(y=C_freq+T_freq, yend=C_freq+T_freq, xend=gencoor), colour= "green") +
  geom_text(aes(y=C_conv,label=paste(action, '\n this happed on ', gencoor)),size=2.5,hjust=-.01, vjust=-.01, angle = 35) +
  geom_point(aes(y=C_conv)) +
  geom_hline(y=0,size=1,color='purple') +
  theme(axis.text.y = element_blank()) + ylab('') + xlab('') +
  labs(title = "Timeline for when what happened")
plot


write.csv(BID_BACS_filtered_AUC_manual_18s, "BID_BACS_filtered_AUC_manual_18s.csv")



#Input file
#Input file
#Input file
#Input file

#Frequencies of nucleotides for treat file
A_freqi <- c()
C_freqi <- c()
G_freqi <- c()
T_freqi <- c()
C_convi <- c()
bkg_C_convi <- c()
bkg_G_convi <- c()
bkg_T_convi <- c()
bkg_A_convi <- c()
for (row in 1:nrow(BACSinput)) {
  cov1 = BACSinput$A[row] + BACSinput$G[row] + BACSinput$C[row] + BACSinput$T[row] 
  A_ratioi <- round(BACSinput$A[row] / cov1 *100, digits = 1)
  A_freqi <- append(A_freqi, A_ratioi)
  C_ratioi <- round(BACSinput$C[row] / cov1*100, digits = 1)
  C_freqi <- append(C_freqi, C_ratioi)
  G_ratioi <- round(BACSinput$G[row] / cov1*100, digits = 1)
  G_freqi <- append(G_freqi, G_ratioi)
  T_ratioi <- round(BACSinput$T[row] / cov1*100, digits = 1)
  T_freqi <- append(T_freqi, T_ratioi)
}
for (row in 1:nrow(BACSinput)) {
  cov3 = BACSinput$C[row] + BACSinput$T[row] 
  A_ratioi <- round(BACSinput$A[row] / cov3 *100, digits = 1)
  bkg_A_convi <- append(bkg_A_convi, A_ratioi)
  C_ratioi <- round(BACSinput$C[row] / cov3 *100, digits = 1)
  bkg_C_convi <- append(bkg_C_convi, C_ratioi)
  G_ratioi <- round(BACSinput$G[row] / cov3 *100, digits = 1)
  bkg_G_convi <- append(bkg_G_convi, G_ratioi)
  T_ratioi <- round(BACSinput$T[row] / cov3 *100, digits = 1)
  bkg_T_convi <- append(bkg_T_convi, T_ratioi)
}
for (row in 1:nrow(BACSinput)) {
  cov2 = BACSinput$C[row] + BACSinput$T[row] 
  C_convrati <- round(BACSinput$C[row] / cov2*100, digits = 1)
  C_convi <- append(C_convi, C_convrati)
}
#Add columns with frequencies
BACSinput$A_freqi = A_freqi
BACSinput$C_freqi = C_freqi
BACSinput$G_freqi = G_freqi
BACSinput$T_freqi = T_freqi
BACSinput$C_convi = C_convi
BACSinput$bkg_A_convi = bkg_A_convi
BACSinput$bkg_C_convi = bkg_C_convi
BACSinput$bkg_G_convi = bkg_G_convi
BACSinput$bkg_T_convi = bkg_T_convi

BACSinput$common <- paste0(BACSinput$gene, "_", BACSinput$gencoor)
#Merge treat+ input
BACStreat_input = merge(x=BACSinput, y=BACStreat, by="common", all.x = T, all.y = T )
write.csv(BACStreat_taoka,"~/BACStreat_input.csv", row.names = FALSE)
BACStreat_input_taoka = merge(x=taoka, y=BACStreat_input, by="common", all.x = T, all.y = T )



ggplot(BACStreat_input_taoka %>% filter(V3 == "Y", gene.x == "18s"), aes(x = gencoor.x)) +
  # Add segments for C_conv and Conv_Taoka
#  geom_segment(aes(xend = gencoor.x, y = 0, yend = Conv_Taoka), color = "darkviolet") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_conv), color = "brown2") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_convi), color = "grey") +
  # Add points for Conv_Taoka
#  geom_point(aes(y = Conv_Taoka), color = "darkviolet", size = 3) +
  # Add points for C_conv
  geom_point(aes(y = C_conv), color = "brown2", size = 3) +
  # Add points for C_convi
  geom_point(aes(y = C_convi), color = "grey", size = 3) +
  # Add text labels for gencoor.x (optional for clarity)
  #geom_text(aes(y = C_conv, label = paste0(gencoor.x)), 
  #         vjust = -0.5, color = "black", size = 2) +
  #geom_text(aes(y = Conv_Taoka, label = paste0(gencoor.x)), 
  #         vjust = -0.5, color = "black", size = 2) +
  # Adjust theme and labels
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position in 18s") +
  ylab("Conversion rate (%)") +
  ggtitle("BACS. Comparison of C_conv and Conv_Taoka Rates")


ggplot(BACStreat_input_taoka %>% filter(V3 == "T", gene.x == "5.8s"), aes(x = gencoor.x, shape = V3 == "Y")) +
  # Add segments for C_conv and Conv_Taoka
  geom_segment(aes(xend = gencoor.x, y = 0, yend = Conv_Taoka), color = "#C3D7A4" ) +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_conv), color = "#FC4E07") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_convi), color = "#4E84C4") +
  geom_point(aes(y = Conv_Taoka), color = "#C3D7A4", size = 3) +
  geom_point(aes(y = C_conv), color = "#FC4E07", size = 3) +
  geom_point(aes(y = C_convi), color = "#4E84C4", size = 3) +
  # Define custom shapes for logical values
  scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 16)) + # Square for FALSE, Circle for TRUE
  
  # Adjust theme and labels
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position in 5.8s") +
  ylab("Conversion rate (%)") +
  ggtitle("BACS: Comparison of C conversion Rates. 5.8s")

#Best and last

BACStreat_input_taoka %>%
  filter(refSeq.y == "T", gene.x == "18s") %>%
  ggplot(aes(x = gencoor.x)) +
  # Add segments with fixed colors
  geom_segment(aes(xend = gencoor.x, y = 0, yend = Conv_Taoka, color = "Conv_Taoka")) +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_conv, color = "C_conv")) +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_convi, color = "C_convi")) +
  
  # Conditional coloring for Conv_Taoka points
  geom_point(aes(y = Conv_Taoka, color = ifelse(V3 == "Y", "V3 == Y", "V3 != Y")), size = 3) +
  
  # Fixed colors for C_conv and C_convi points
  geom_point(aes(y = C_conv, color = "C_conv"), size = 3) +
  geom_point(aes(y = C_convi, color = "C_convi"), size = 3) +
  
  # Define color mappings and legend labels
  scale_color_manual(
    values = c(
      "Conv_Taoka" = "#C3D7A4",
      "C_conv" = "#FC4E07",
      "C_convi" = "#4E84C4",
      "V3 == Y" = "#C3D7A4",
      "V3 != Y" = "#009E73"
    ),
    name = "Legend",
    breaks = c("C_conv", "C_convi", "V3 == Y", "V3 != Y"),
    labels = c("C_conv_treat", "C_conv_input", "Taoka, Y", "Taoka, non-Y")
  ) +
  
  # Minimal theme and labels
  theme_minimal() +
  labs(
    title = "BACS.C Conversion. 18s",
    x = "Genomic Coordinate",
    y = "Values"
  )


ggplot(BIDtreat_rRNA_Del_taoka %>% filter(refSeq.y == "T", gene.x == "18s"), 
       aes(x = gencoor.x, shape = V3 == "Y")) +
  # Add segments for C_conv and Conv_Taoka with color mapping
  geom_segment(aes(xend = gencoor.x, y = 0, yend = Conv_Taoka, color = "Conv_Taoka")) +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = Del_conv, color = "Del_conv")) +
  
  # Points for Conv_Taoka and Del_conv with color mapping
  geom_point(aes(y = Conv_Taoka, color = "Conv_Taoka"), size = 3) +
  geom_point(aes(y = Del_conv, color = "Del_conv"), size = 3) +
  
  # Define custom colors for Conv_Taoka and Del_conv
  scale_color_manual(
    values = c("Conv_Taoka" = "#C3D7A4", "Del_conv" = "#FC4E07"),
    name = "Conversion Type"  # Title for the legend
  ) +
  
  # Define custom shapes for logical values
  scale_shape_manual(
    values = c("FALSE" = 15, "TRUE" = 16),
    name = "Deletion Present",
    labels = c("No", "Yes")
  ) +
  
  # Adjust theme and labels
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Position in 18s") +
  ylab("Conversion rate (%)") +
  ggtitle("BID: Comparison of Del_conv and Conv_Taoka Rates")






















library(ggplot2)
library(dplyr)

BACStreat_input_taoka %>%
  filter(refSeq.y == "T", gene.x == "18s") %>%
  ggplot(aes(x = gencoor.x, color = V3)) +  # Map 'shape' here
  geom_segment(aes(xend = gencoor.x, y = 0, yend = Conv_Taoka), color = "#C3D7A4") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_conv), color = "#FC4E07") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_convi), color = "#4E84C4") +
  geom_point(aes(y = Conv_Taoka), size = 3, color = "#C3D7A4") +
  geom_point(aes(y = C_conv), size = 3, color = "#FC4E07") +
  geom_point(aes(y = C_convi), size = 3, color = "#4E84C4") +
  scale_shape_manual(values = c(`FALSE` = 15, `TRUE` = 16)) +  # Correct manual shape assignment
  theme_minimal() +
  labs(
    title = "Comparison of Conv_Taoka, C_conv, and C_convi",
    x = "Genomic Coordinate",
    y = "Values"
  )

BACStreat_input_taoka %>%
  filter(V3 == "Y", gene.x == "28s") %>%
  ggplot(aes(x = gencoor.x)) +  # Map 'shape' here
  geom_segment(aes(xend = gencoor.x, y = 0, yend = Conv_Taoka), color = "#C3D7A4") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_conv), color = "#FC4E07") +
  geom_segment(aes(xend = gencoor.x, y = 0, yend = C_convi), color = "#4E84C4") +
  geom_point(aes(y = Conv_Taoka), size = 3, color = "#C3D7A4") +
  geom_point(aes(y = C_conv), size = 3, color = "#FC4E07") +
  geom_point(aes(y = C_convi), size = 3, color = "#4E84C4") +
  scale_shape_manual(values = c(`FALSE` = 15, `TRUE` = 16)) +  # Correct manual shape assignment
  theme_minimal() +
  labs(
    title = "BACS. C conversion. 28s",
    x = "Genomic Coordinate",
    y = "Values"
  )



colnames(treat_taoka)

library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)



# Plot the data
ggplot(treat_T_taoka%>%filter(V3== "Y"), aes(x = gencoor.y, y = bkg_C_conv, color = V3)) +
  geom_line() +
  theme_light() +
  xlab("Genomic Coordinates") +
  ylab("Background Conversion Rate") +
  ggtitle("Background Conversion Rates for T Reference Sequence")

#Background plot 
ggplot(BACStreat_input_taoka%>% filter(refSeq.x == "T", gene.x == "5.8s"), aes(x = gencoor.x, y = bkg_C_conv, color = V3)) +
  geom_point(size = 3) +
  theme_light() +
  xlab("Genomic Coordinates") +
  ylab("Conversion Rate (%)") +
  ggtitle("BACS. Conversion Rates for T Reference Sequence. 5.8s")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  ) +
  scale_color_manual(
    values = c(
      "m1acp3Y" = "purple", 
      "pU" = "blue", 
      "Um" = "green", 
      "Y" = "red", 
      "NA" = "black"
    )
  )+
  ylim(0, 100)

ggplot(treat_18s_taoka%>% filter(refSeq.y == "T"), aes(x = gencoor, y = bkg_A_conv, color = V3)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("Genomic Coordinates") +
  ylab("Conversion Rate (%)") +
  ggtitle("Conversion Rates for T Reference Sequence")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
  ) +
  scale_color_manual(
    values = c(
      "m1acp3Y" = "red", 
      "pU" = "blue", 
      "Um" = "green", 
      "Y" = "purple", 
      "NA" = "black"
    )
  )+
  ylim(0, 100)


ggplot(treat_18s%>% filter(refSeq == "T"), aes(x =gencoor, y = bkg_T_conv, color = variable)) +
  geom_line()

plot <- ggplot(treat_18s_T, aes(x=gencoor)) +
  #Three geom_segments 
  geom_segment(aes(y=C_conv,yend=C_conv,xend=gencoor), colour= "blue") +
  geom_segment(aes(y=C_freq,yend=C_freq,xend=gencoor), colour= "red") +
  geom_segment(aes(y=C_freq+T_freq, yend=C_freq+T_freq, xend=gencoor), colour= "green") +
  geom_text(aes(y=C_conv,label=paste(action, '\n this happed on ', gencoor)),size=2.5,hjust=-.01, vjust=-.01, angle = 35) +
  geom_point(aes(y=C_conv)) +
  geom_hline(y=0,size=1,color='purple') +
  theme(axis.text.y = element_blank()) + ylab('') + xlab('') +
  labs(title = "Timeline for when what happened")
plot


write.csv(BID_BACS_filtered_AUC_manual_18s, "BID_BACS_filtered_AUC_manual_18s.csv")


library(ggplot2)
library(dplyr)


