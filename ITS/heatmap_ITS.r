library("tidyverse")
library ("phyloseq")
library("ComplexHeatmap")
library("circlize")


dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq_ITS.RData")



############################
############################
# Metagenome 
############################
############################

##############################################
# Merge samples (Organ_material)
##############################################

physeq = physeq_ITS
physeq = merge_samples(physeq, "Organ_material")
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)


physeq = psmelt(physeq) %>%
	group_by(Genus, Sample) %>%
	summarise(Abundance = sum(Abundance)) 


Heatmap = pivot_wider(physeq, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
Heatmap = as.data.frame(Heatmap)
row.names(Heatmap) = Heatmap[,1]
Heatmap = Heatmap[,-1]

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 10) %>% # Top 10 taxa of genus level in heatmap
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Figure 4. Heatmap Organ_material (ITS).png",  width = 20, height = 15, units = "in", res = 300)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "",
	column_title_gp = gpar(fontsize = 20),
	cluster_columns = FALSE,
	column_order = c("Branch","Needles","Shavings"),
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_gp = gpar(fontsize = 20),
	column_names_rot = 90,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 20), 
		title_gp = gpar(fontsize = 20)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 17))
})

draw(Hmap)

dev.off()



##############################################
# Merge samples (Season)
##############################################

physeq = physeq_ITS
physeq = merge_samples(physeq, "Season")
physeq = transform_sample_counts(physeq, function(x) x / sum(x) * 100)


physeq = psmelt(physeq) %>%
	group_by(Genus, Sample) %>%
	summarise(Abundance = sum(Abundance)) 


Heatmap = pivot_wider(physeq, id_cols = "Genus", names_from = "Sample", values_from = "Abundance")
Heatmap = as.data.frame(Heatmap)
row.names(Heatmap) = Heatmap[,1]
Heatmap = Heatmap[,-1]

Heatmap_top = Heatmap
Heatmap_top[Heatmap_top[] == 0] <- NA

for(loc in colnames(Heatmap_top)){
Heatmap_top[,loc] = rank(-Heatmap_top[,loc], ties.method = "first", na.last = "keep")
}
Heatmap_top[is.na(Heatmap_top)] <- 10000

Heatmap_top = Heatmap_top %>%
	mutate(Min = do.call(pmin, select_if(., is.numeric)))

Heatmap_top = Heatmap_top %>%
	filter(Min <= 10) %>% # Top 10 taxa of genus level in heatmap
	select(!Min)

Heatmap = filter(Heatmap, rownames(Heatmap) %in% rownames(Heatmap_top))

Heatmap_percent = as.matrix(Heatmap)
Heatmap_percent[Heatmap_percent[] == 0] <- 1000
Heatmap_percent[Heatmap_percent[] < 0.1] <- 10000
Heatmap_percent[] = scales::percent(Heatmap_percent, scale = 1, accuracy = 0.01)
Heatmap_percent[Heatmap_percent[] == scales::percent(1000, scale = 1, accuracy = 0.01)] <- NA
Heatmap_percent[Heatmap_percent[] == scales::percent(10000, scale = 1, accuracy = 0.01)] <- "<0.10%"

png("./Figure 7b. Heatmap Seasons (ITS).png",  width = 15, height = 9, units = "in", res = 500)

col_fun = colorRamp2(breaks = c(0, max(Heatmap)), colors = c("#ffffff", "#ff0000"))

Hmap = Heatmap(as.matrix(Heatmap), 
	name = "Realative abundance, %", 
	col = col_fun,
	column_title = "",
	column_title_gp = gpar(fontsize = 20),
	cluster_columns = FALSE,
	column_order = c("Summer","Autumn","Winter"),
	row_names_max_width = max_text_width("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
	row_names_gp = gpar(fontsize = 20, fontface = "italic"),
	column_names_gp = gpar(fontsize = 20),
	column_names_rot = 90,
	column_names_centered = TRUE,
	heatmap_legend_param = list(legend_height = unit(5, "cm"), 
		grid_width = unit(1, "cm"), 
		labels_gp = gpar(fontsize = 20), 
		title_gp = gpar(fontsize = 20)),
	cell_fun = function(j, i, x, y, w, h, col) {
grid.text(Heatmap_percent[i, j], x, y, gp = gpar(fontsize = 17))
})

draw(Hmap)
dev.off()
