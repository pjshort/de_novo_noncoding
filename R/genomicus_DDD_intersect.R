# match the ddd de novos with target genes using genomicus annotation

de_novos = read.table("../data/de_novo_filtered.txt", sep = "\t", header = TRUE)
genomicus_CREs = read.table("../data/CREs_intersect_REG.list", header = TRUE, sep = "\t")

