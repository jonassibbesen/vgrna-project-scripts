
# impute_hla_cds_alleles.R

rm(list=ls())

library(hlaseqlib)

args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args) != 3) {
  
  print("Usage: Rscript impute_hla_cds_alleles.R <gene_name> <imgthla_db_folder> <output_name>")
  quit(status = 1)
}

hla_alleles <- hla_compile_index(args[1], args[2])

print(nrow(hla_alleles))

write.table(hla_alleles, args[3], row.name = F, quote = F, sep = "\t")


