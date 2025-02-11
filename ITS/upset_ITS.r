library ("RColorBrewer")
library("tidyverse")
library ("phyloseq")
library("circlize")
library("ComplexHeatmap")



dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq_ITS.RData")



############################
############################
# Metagenome 
############################
############################

##############################################
# Filter low abundand taxa (Filter above species level taxa with less than 0.5% relative abundance among all samples)
##############################################

filter_genera = transform_sample_counts(physeq_ITS, function(x) x / sum(x) * 100)

filter_genera = psmelt(filter_genera) %>%
	group_by(Genus, Sample) %>%
	summarise(Abundance = sum(Abundance)) %>%
	group_by(Genus) %>%
	mutate(Genus = ifelse(max(Abundance) > 0.5, Genus, "Other")) %>% 
	select(Genus) %>%
	filter(Genus != "Other") %>%
	unique()

##############################################
# Merge samples (Organ_material)
##############################################

physeq = merge_samples(physeq_ITS, "Organ_material")
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)

physeq_main = psmelt(physeq) %>%
	group_by(Genus, Sample) %>%
	summarise(Abundance = sum(Abundance)) %>%
	group_by(Genus) %>%
	mutate(Genus = ifelse(Genus %in% filter_genera$Genus, Genus, "Other")) %>% 
	filter(Genus != "Other") %>%
	select(Genus, Sample, Abundance)



pre_upset_data = pivot_wider(physeq_main, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
pre_upset_data = as.data.frame(pre_upset_data)
row.names(pre_upset_data) = pre_upset_data[,1]
pre_upset_data = pre_upset_data[,-1]

upset_data = pre_upset_data

head(upset_data)
upset_data[] = ifelse(upset_data[] > 0, 1, 0) 
upset_data = upset_data[!apply(upset_data == 0, 1, all),]


upset_data_for_suppl = upset_data

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)
for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}
upset_data_for_suppl = t(upset_data_for_suppl)[-1,]
upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}
Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Genus", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Table S7.txt")
Upset_genus_.txt = left_join(Upset_genus_.txt, tax_suppl, by = "Genus")
Upset_genus_.txt = Upset_genus_.txt[,c("Where present", "Amplicon sequence variant", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")]
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Table S6.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 3c. Upset Organ material (ITS)
png("Figure 3c. Upset_Organ_material (ITS).png",  width = 8, height = 6, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 22, fontface="italic"),
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 12)),
		annotation_name_gp = gpar(fontsize = 12),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 16)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 12)),
		annotation_name_gp = gpar(fontsize = 12),
		numbers_gp = gpar(fontsize = 16)))
draw(upset_plot)
dev.off()


##############################################
# Merge samples (Season)
##############################################


physeq = merge_samples(physeq_ITS, "Season")
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)

physeq_main = psmelt(physeq) %>%
	group_by(Genus, Sample) %>%
	summarise(Abundance = sum(Abundance)) %>%
	group_by(Genus) %>%
	mutate(Genus = ifelse(Genus %in% filter_genera$Genus, Genus, "Other")) %>% 
	filter(Genus != "Other") %>%
	select(Genus, Sample, Abundance)



pre_upset_data = pivot_wider(physeq_main, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
pre_upset_data = as.data.frame(pre_upset_data)
row.names(pre_upset_data) = pre_upset_data[,1]
pre_upset_data = pre_upset_data[,-1]

upset_data = pre_upset_data

head(upset_data)
upset_data[] = ifelse(upset_data[] > 0, 1, 0) 
upset_data = upset_data[!apply(upset_data == 0, 1, all),]


upset_data_for_suppl = upset_data

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)
for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}
upset_data_for_suppl = t(upset_data_for_suppl)[-1,]
upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}
Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Genus", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Table S12.txt")
Upset_genus_.txt = left_join(Upset_genus_.txt, tax_suppl, by = "Genus")
Upset_genus_.txt = Upset_genus_.txt[,c("Where present", "Amplicon sequence variant", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")]
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Table S11.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 6d. Upset Season (ITS)
png("Figure 6d. Upset_Season (ITS).png",  width = 8, height = 6, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 22, fontface="italic"),
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 12)),
		annotation_name_gp = gpar(fontsize = 12),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 16)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 12)),
		annotation_name_gp = gpar(fontsize = 12),
		numbers_gp = gpar(fontsize = 16)))
draw(upset_plot)
dev.off()




##############################################
# Merge samples (Plant)
##############################################


physeq = merge_samples(physeq_ITS, "Plant")
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)

physeq_main = psmelt(physeq) %>%
	group_by(Genus, Sample) %>%
	summarise(Abundance = sum(Abundance)) %>%
	group_by(Genus) %>%
	mutate(Genus = ifelse(Genus %in% filter_genera$Genus, Genus, "Other")) %>% 
	filter(Genus != "Other") %>%
	select(Genus, Sample, Abundance)



pre_upset_data = pivot_wider(physeq_main, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
pre_upset_data = as.data.frame(pre_upset_data)
row.names(pre_upset_data) = pre_upset_data[,1]
pre_upset_data = pre_upset_data[,-1]

upset_data = pre_upset_data

head(upset_data)
upset_data[] = ifelse(upset_data[] > 0, 1, 0) 
upset_data = upset_data[!apply(upset_data == 0, 1, all),]


upset_data_for_suppl = upset_data

upset_data_for_suppl = t(upset_data_for_suppl)
upset_data_for_suppl = cbind(rownames(upset_data_for_suppl), upset_data_for_suppl)
for(y in c(2:length(colnames(upset_data_for_suppl)))) {
	for(x in c(1:length(rownames(upset_data_for_suppl)))) {
		if(upset_data_for_suppl[x,y] == 1) {upset_data_for_suppl[x,y] = upset_data_for_suppl[x,1]} 
			else {upset_data_for_suppl[x,y] = 100}
	}
}
upset_data_for_suppl = t(upset_data_for_suppl)[-1,]
upset_data_for_suppl_rows = rownames(upset_data_for_suppl)
Upset_genus_.txt = NULL
for(x in c(1:length(rownames(upset_data_for_suppl)))) {
	pasted_cols = upset_data_for_suppl[x, 1:length(colnames(upset_data_for_suppl))]
	pasted_cols = pasted_cols[!pasted_cols == 100]
	pasted_cols = paste(pasted_cols, collapse = ", ")
	Upset_genus_.txt = rbind(Upset_genus_.txt,pasted_cols)
}
Upset_genus_.txt = cbind(upset_data_for_suppl_rows, Upset_genus_.txt)
colnames(Upset_genus_.txt) = c("Genus", "Where present")
Upset_genus_.txt = as.data.frame(Upset_genus_.txt)
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Table S16.txt")
Upset_genus_.txt = left_join(Upset_genus_.txt, tax_suppl, by = "Genus")
Upset_genus_.txt = Upset_genus_.txt[,c("Where present", "Amplicon sequence variant", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")]
write.table(Upset_genus_.txt, sep = "\t", row.names = FALSE, quote = FALSE, file="Table S15.txt")
#

upset_data = make_comb_mat(upset_data)

upset_color_palette = brewer.pal(8, "Set1")
nb.cols = length(comb_name(upset_data))
upset_colors = colorRampPalette(upset_color_palette)(nb.cols)

# Figure 9d. Upset Plant (ITS)
png("Figure 9d. Upset_Plant(ITS).png",  width = 10, height = 8, units = "in", res = 300)
upset_plot = UpSet(upset_data, 
	comb_col = upset_colors, 
    row_names_gp = gpar(fontsize = 22, fontface="italic"),
	top_annotation = upset_top_annotation(upset_data, 
		numbers_rot = 0, 
		axis_param = list(gp = gpar(fontsize = 12)),
		annotation_name_gp = gpar(fontsize = 12),
		add_numbers = TRUE, 
		numbers_gp = gpar(fontsize = 16)),
    right_annotation = upset_right_annotation(upset_data, 
		add_numbers = TRUE, 
		axis_param = list(gp = gpar(fontsize = 12)),
		annotation_name_gp = gpar(fontsize = 12),
		numbers_gp = gpar(fontsize = 16)))
draw(upset_plot)
dev.off()




