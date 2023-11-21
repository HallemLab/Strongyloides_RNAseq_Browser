# Script for saving 

species_list <- tibble(species = c('Ss', 'Sr', 'Sv', 'Sp'))

for (x in species_list$species) {
# Import a variance-stabilized DEGList created by voom transformation command.
# Outputs: E = normalized CPMexpression values on the log2 scale
load(file = paste0("../Data/",x,"_vDGEList"))
v.DEGList.filtered.norm <- v.DEGList.filtered.norm

target.contrast.options <- v.DEGList.filtered.norm$targets$group

# Import a tidy dataframe containing gene annotations for all genes in the genome (including those that are excluded from this database.)
load(file = paste0("../Data/",x,"_geneAnnotations"))
annotations <- as_tibble(annotations, rownames = "geneID")

# Parse vDGEList into a tibble containing Log2CPM information
Log2CPM<-v.DEGList.filtered.norm$E %>%
    as_tibble(rownames = "geneID")%>%
    setNames(nm = c("geneID", 
                    as.character(v.DEGList.filtered.norm$targets$group))) %>%
    pivot_longer(cols = -geneID,
                 names_to = "life_stage", 
                 values_to = "log2CPM") %>%
    group_by(geneID, life_stage)

diffGenes.df <- v.DEGList.filtered.norm$E %>%
    as_tibble(rownames = "geneID", .name_repair = "unique")

## Fit a linear model to the data
fit <- lmFit(v.DEGList.filtered.norm, v.DEGList.filtered.norm$design)
save(fit, file = paste0("../Data/", x, "_lmFit"))
}

