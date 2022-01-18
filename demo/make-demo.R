library(phyloseq)

data("GlobalPatterns")

food <-
  import_biom(
    BIOMfilename = "demo/chaillou.biom",
    treefilename = "demo/chaillou.nwk",
    parseFunction = parse_taxonomy_greengenes
  )
sample_data(food) <-
  read.table("demo/chaillou.tsv",
             header = TRUE,
             row.names = 1,
             sep = "\t")

# kinetic <-
#   import_biom(BIOMfilename = "mach.biom", treefilename = "mach.nwk")
# sample_data(kinetic) <-
#   read.table("mach.tsv",
#              header = TRUE,
#              row.names = 1,
#              sep = "\t")
# colnames(tax_table(kinetic)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# tax_table(kinetic)[grep("unknown ", tax_table(kinetic))] <- NA
# 
# ravel <-
#   import_biom(BIOMfilename = "ravel.biom")
# sample_data(ravel) <-
#   read.table("ravel.tsv",
#              header = TRUE,
#              row.names = 1,
#              sep = "\t")
# colnames(tax_table(ravel)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# tax_table(ravel)[grep("unknown ", tax_table(ravel))] <- NA
# 
# soil <-
#   import_biom(BIOMfilename = "morton.biom", treefilename = "morton.tree")
# sample_data(soil) <-
#   read.table("morton.csv",
#              header = TRUE,
#              row.names = 1,
#              sep = "\t")
# colnames(tax_table(soil)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# tax_table(soil)[grep("unknown ", tax_table(soil))] <- NA
# 
# biorare <-
#   import_biom(BIOMfilename = "biorare.biom", treefilename = "biorare.nwk")
# sample_data(biorare) <-
#   read.table("biorare.tsv", header = TRUE, row.names = 1)
# colnames(tax_table(biorare)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# tax_table(biorare)[grep("unknown ", tax_table(biorare))] <- NA
# 
# cellulose <- import_biom(BIOMfilename = "cellulose.biom")
# sample_data(cellulose) <- read.table("cellulose.tsv", header = TRUE, row.names = 1)

# save(GlobalPatterns, food, kinetic, soil, ravel, biorare, file = "demo.RData")

saveRDS(food, file = "demo/food.RDS")
saveRDS(GlobalPatterns, file = "demo/GlobalPatterns.RDS")