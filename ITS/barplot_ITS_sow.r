library("tidyverse")
library ("phyloseq")

colors_16s = c("#69d2e7", "#a7dbd8", "#e0e4cc", "#f38630", "#fa6900", "#fc9d9a", "#f9cdad", 
"#c8c8a9", "#83af9b", "#ecd078", "#d95b43", "#542437", "#53777a", 
"#4ecdc4", "#c7f464", "#ff6b6b", "#c44d58", "#774f38", "#e08e79", "#f1d4af", 
"#c5e0dc", "#018BFF", "#01FF61", "#FFD101", 
"#C02727", "#C177C7", "#EBBAE3", "#BAEBC4", "#397DD7", "#6FE8D9")

colors_ITS = c("#018BFF", "#01FDFF", "#01FF61", "#F9FF01", "#FFD101", "#C08227", 
"#C02727", "#ecd078", "#C177C7", "#EBBAE3", "#BAEBC4", "#397DD7", "#6FE8D9", 
"#69d2e7", "#a7dbd8", "#e0e4cc", "#f38630", "#fa6900", "#fe4365", "#fc9d9a", "#f9cdad", 
"#c8c8a9", "#83af9b", "#65216B", "#d95b43", "#c02942", "#542437", "#53777a", 
"#4ecdc4", "#c7f464", "#ff6b6b", "#c44d58", "#774f38", "#e08e79", "#f1d4af", 
"#c5e0dc")


mycolors = colors_ITS

dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq_ITS.RData")

# For classes ordering and coloring in bar plots
physeq_factors = psmelt(physeq_ITS) %>%
	select(Class, Sample, Abundance) %>%
	group_by(Class, Sample) %>%
	summarise(Abundance = sum(Abundance))

physeq_factors = pivot_wider(physeq_factors, id_cols = "Class", names_from = "Sample", values_from = "Abundance")
physeq_factors = as.data.frame(physeq_factors)
row.names(physeq_factors) = physeq_factors[,1]
physeq_factors = physeq_factors[,-1]

physeq_factors$Mean = rowMeans(physeq_factors[])
physeq_factors$Class = row.names(physeq_factors)
physeq_factors = physeq_factors[,c("Class", "Mean")]
physeq_factors = pivot_longer(physeq_factors, cols = !Class, names_to = "Sample", values_to = "Abundance")
physeq_factors = physeq_factors[order(physeq_factors$Abundance),]
physeq_factors = physeq_factors$Class
# physeq_factors = physeq_factors[!physeq_factors == "Other"]
physeq_factors = cbind(physeq_factors, c("#aed7ea" ,rev(mycolors[1:length(physeq_factors)-1])))
colnames(physeq_factors) = c("Class", "color")



#


############################
############################
# Culture-dependent method
############################
############################

##############################################
# Merge samples (Organ) and add "Sum" column 
##############################################


MpsplS = physeq_ITS_sow
sample_data(MpsplS)$Organ_material = "Average_mean"
MpsplS = merge_samples(MpsplS, "Organ_material")

psPl = physeq_ITS_sow
psPl = merge_samples(psPl, "Organ_material")

Spspl = merge_phyloseq(MpsplS, psPl)
#

# For read counts above bars in bar plots
physeq_read_counts = psmelt(psPl) %>%
	select(Class, Sample, Abundance) %>%
	group_by(., Class, Sample) %>%
	summarise(., Abundance = sum(Abundance))
For_barplot = pivot_wider(physeq_read_counts, id_cols = "Class", names_from = "Sample", values_from = "Abundance")
For_barplot = For_barplot[,-1]
For_barplot = colSums(For_barplot)
For_barplot = t(as.matrix(For_barplot))
row.names(For_barplot) = "Sum"
For_barplot = as.data.frame(For_barplot)
For_barplot = pivot_longer(For_barplot, cols = c(colnames(For_barplot)), names_to = "Sample", values_to = "Abundance")
colnames(For_barplot) = c("Sample", "Sum")
For_barplot[,2] = apply(For_barplot[,2], 1, function(x) formatC(x, big.mark=",", format = "d"))
#


Spspl = transform_sample_counts(Spspl, function(x) x / sum(x) * 100)


Spspl = psmelt(Spspl) %>%
 	group_by(Class, Sample) %>%
	summarise(Abundance = sum(Abundance))


Barplot_factors = physeq_factors[is.element(physeq_factors[,"Class"], unique(Spspl$Class)),]
Spspl$Class = factor(Spspl$Class, levels = Barplot_factors[, "Class"])

Spspl = left_join(Spspl, For_barplot, by="Sample")

Spspl$Sample = factor(Spspl$Sample, levels = c("Branch","Needles","Shavings","Average_mean"))

# Figure. Class level barplot among plant organs (ITS)
ggplot(Spspl) +
 	geom_col(aes(x = Sample, y = Abundance, fill = Class)) +
	scale_fill_manual(values = Barplot_factors[, "color"]) +
	guides(fill = guide_legend(ncol = 1, byrow = TRUE)) + 
	geom_text(aes(x = Sample, y = Abundance, 
		label = ifelse(Abundance > 2, scales::percent(Abundance, scale = 1, accuracy = 1), ""), 
		group = Class), 
		inherit.aes = FALSE, position = position_stack(vjust = 0.5), 
		size = 11, col = "black") +
	geom_text(aes(x = Sample, y = 102, label=Sum), size = 9) +
	theme(legend.text = element_text(size = 27, colour = "black"),
		legend.title = element_text(size = 25),
		axis.text.x = element_text(size=27, angle=0, vjust = 0.5),
		axis.text.y = element_text(size=25),
		axis.title.y = element_text(size=22)) +
	labs(x = "", y = "Relative abundance, %")

ggsave("Figure 3b. barplot_Class_Organ_ITS (sow).png", width = 17, height = 15)

#

