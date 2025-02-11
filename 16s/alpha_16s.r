library("tidyverse")
library ("phyloseq")
library("vegan")

dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./16s/physeq_16s.RData")


alpha_diversity = rarefy_even_depth(physeq_16s)
alpha_diversity = estimate_richness(alpha_diversity, split = TRUE, measures = NULL)
alpha_diversity$SRA = rownames(alpha_diversity)
alpha_diversity = left_join(sampledata_16s, alpha_diversity, by = "SRA")


# Figure 5a. Shannon′s alpha diversity among organs (16s)
alpha_diversity$Organ_material = factor(alpha_diversity$Organ_material, level = c("Branch","Needles","Shavings"))
ggplot(data = alpha_diversity, aes(x = Organ_material, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") + 
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("./Figure 5a. Shannon′s alpha diversity among organs (16s).png", width = 7, height = 10)
#


# Pairwise Wilcoxon rank sum test (organs) (16s)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Organ_material, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test (organs) (16s)")
Shannon$p.value
write.table(Shannon$p.value, sep = "\t", quote = FALSE, file="Pairwise Wilcoxon rank sum test (organs) (16s).txt")



# Figure 8a. Shannon′s alpha diversity among seasons (16s)
alpha_diversity$Season = factor(alpha_diversity$Season, level = c("Summer","Autumn","Winter"))
ggplot(data = alpha_diversity, aes(x = Season, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("./Figure 8a. Shannon′s alpha diversity among seasons (16s).png", width = 7, height = 10)
#
dev.off()

# Pairwise Wilcoxon rank sum test (seasons) (16s)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Season, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test (seasons) (16s)")
Shannon$p.value
write.table(Shannon$p.value, sep = "\t", quote = FALSE, file="Pairwise Wilcoxon rank sum test (seasons) (16s).txt")


# Figure 10a. Shannon′s alpha diversity among locations (16s)
alpha_diversity$location = factor(alpha_diversity$location, level = c("Primorsky_Krai","Sakhalin_Oblast"))
ggplot(data = alpha_diversity, aes(x = location, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("./Figure 10a. Shannon′s alpha diversity among locations (16s).png", width = 7, height = 10)
#
dev.off()

# Pairwise Wilcoxon rank sum test (locations) (16s)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$location, p.adjust.method = 'none')
print("Pairwise Wilcoxon rank sum test (locations) (16s)")
Shannon$p.value
write.table(Shannon$p.value, sep = "\t", quote = FALSE, file="Pairwise Wilcoxon rank sum test (locations) (16s).txt")



