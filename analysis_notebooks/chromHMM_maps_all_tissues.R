# test the MAPS difference between quiescent and all other states in each tissue type

source("~/software/dddMAPS/dddMAPS/MAPS.R")
load("~/software/dddMAPS/data/DDD_4k_parents_synonymous_maps_lm.RData")
library(stringr)

mu_snp <- read.table("~/reference_data/forSanger_1KG_mutation_rate_table.txt", header=TRUE)
gencode = read.table("~/reference_data/gencode.v19.CDS.probe_overlap.min10_coverage.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
noncoding_intervals = read.table("~/reference_data/noncoding_control_and_functional.min10_coverage.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)  # only needed for noncoding analysis

sequences = rbind(gencode[,c("chr", "start", "stop", "seq")], noncoding_intervals[,c("chr", "start", "stop", "seq")])

get_MAPS_for_tissue_chromHMM = function(unaff_parent_variants, chromHMM_bed , maps_lm) {
  # load the BED file
  chromHMM_15state = read.table(gzfile(chromHMM_bed), header = FALSE, sep = "\t")
  colnames(chromHMM_15state) = c("chr", "start", "stop", "chromHMM")
  
  unaff_parent_variants$chromHMM_fetal_brain = get_chromHMM(unaff_parent_variants, chromHMM_15state)
  
  print(sprintf("Working on chromHMM MAPS for: %s", chromHMM_bed))
  
  m = maps_adjust(unaff_parent_variants, unaff_parent_variants$chromHMM_fetal_brain, maps_lm = maps_lm)
  counts = table(unaff_parent_variants$chromHMM_fetal_brain)
  counts = counts[names(m$ps_adjusted)]

  ps_quiescent = m$ps_adjusted["15_Quies"]
  se_quiescent = m$standard_error["15_Quies"]
  counts_quiescent = counts["15_Quies"]
  
  tissue = str_match(chromHMM_bed, "E[0-1][0-9][0-9]")[1]
  
  pval_df = data.frame(tissue = str_match(chromHMM_bed, "E[0-1][0-9][0-9]")[1])  # to do - use grepl or sub to get E### for tissue
  maps_points_df = data.frame(tissue = str_match(chromHMM_bed, "E[0-1][0-9][0-9]")[1])  # to do - use grepl or sub to get E### for tissue
  
  for (i in seq_along(names(m$ps_adjusted))) {
    q_count = counts
    test_statistic = (ps_quiescent - m$ps_adjusted[i])/sqrt(se_quiescent^2 + m$standard_error[i]^2)
    pval = -log10(as.numeric(pt(-abs(test_statistic), counts_quiescent + counts[i] - 2)))
    maps_points = (m$ps_adjusted[i] - ps_quiescent)
    if (is.na(pval)) {
      pval = 0
    }
    
    pval_df[,names(m$ps_adjusted)[i]] = pval
    pval_df$score = "pval"
    maps_points_df[,names(m$ps_adjusted)[i]] = maps_points
    maps_points_df$score = "MAPS_points"
    
    df = rbind(pval_df, maps_points_df)
  }
  
  return(df)
}

file_list = list.files(pattern = "/lustre/scratch113/projects/ddd/users/ps14/REP/chromHMM/E[0-1][0-9][0-9]_15_coreMarks_mnemonics.bed.gz")

noncoding_functional_elements = read.table("~/reference_data/noncoding_elements.probe_overlap.min10_coverage.txt", header = TRUE, sep = "\t")
conserved_elements = subset(noncoding_functional_elements, annotation == "Conserved")

unaff_parent_variants = read.table("/lustre/scratch113/projects/ddd/users/ps14/parental_unaffected/unaffected_parent_alleles_all_chromosomes.txt", header = TRUE, sep = "\t")
unaff_parent_variants_conserved = filter_with_bed(unaff_parent_variants, conserved_elements)

all_tissues = lapply(file_list, function(f) get_MAPS_for_tissue_chromHMM(unaff_parent_variants_conserved, f, maps_lm))
df = do.call(rbind, all_tissues)

write.table(df, file = "/lustre/scratch113/projects/ddd/users/ps14/REP/MAPS_all_tissues_chromHMM_15_state.txt", col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)

