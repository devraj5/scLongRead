# SQANTI3 Transcript Filtering Script
# Adapted from Towfique Raj paper: https://doi.org/10.1038/s41588-025-02099-0

library(tidyr)
library(rtracklayer)
library(Biostrings)
library(readr)
library(stringr) 
library(dplyr) 

# CONFIGURATION 
# Setting file paths
sqanti_file <- "/data/project/rizzardilab/drbhatta/long_read/AD_EDA/AD_splice_validation/sqanti3_optimized/optimized_sqanti3_results_classification.txt"
fasta_file <- "/data/project/rizzardilab/drbhatta/long_read/AD_EDA/AD_splice_validation/sqanti3_optimized/optimized_sqanti3_results_corrected.fasta"
gff_file <- "/data/project/rizzardilab/drbhatta/long_read/AD_EDA/AD_splice_validation/sqanti3_optimized/optimized_sqanti3_results_corrected.gtf"


# Filtering parameters
filter_mono_exonic <- TRUE  # Set to TRUE to remove all mono-exonic transcripts

#  LOADING DATA 
message(" * Reading SQANTI classification file")
pre <- read_tsv(sqanti_file)

#  CALCULATING FILTERING METRICS 
# Calculating N_A_downstream_TTS (missing from your file)
message(" * Calculating number of A's in downstream sequence")
pre$N_A_downstream_TTS <- stringr::str_count(pre$seq_A_downstream_TTS, "A")

# Calculating the Roy and Chanfreau score for polyA filtering
# This counts adenines in the downstream 19 bp, counting any adenines in first 6 bp twice
calculate_rc_score <- function(seqs) {
  a_start <- stringr::str_count(stringr::str_sub(seqs, start = 1, end = 6), "A")
  a_total <- stringr::str_count(stringr::str_sub(seqs, start = 1, end = 19), "A")
  a_score <- a_start + a_total
  return(a_score)
}

pre$rc_score <- calculate_rc_score(pre$seq_A_downstream_TTS)


# APPLY FILTERING RULES 
# Based on SQANTI3 default filtering rules:
# 1. FSM: Keep 
# 2. Non-FSM: Keep if:
#    - No RT-switching junction
#    - Has annotated TTS OR low adenosine content at TTS
#    - All junctions canonical or sufficient short read support

message(" * Applying SQANTI filtering rules")
post <- pre %>% 
  mutate(filter_pass = case_when(
    # FSM rule: Keep all FSMs unless they have unreliable 3' ends
    structural_category == "full-splice_match" ~ TRUE,
    
    # Non-FSM rule: Must pass all criteria
    structural_category != "full-splice_match" & 
      RTS_stage == FALSE &  # No RT-switching junction
      (
        diff_to_gene_TTS == 0 |  # Has an annotated TTS
          (N_A_downstream_TTS < 6 & perc_A_downstream_TTS < 60 & rc_score <= 15)  # Or TTS has low A content
      ) ~ TRUE,
    
    TRUE ~ FALSE
  ))

# Filtering mono-exonic transcripts 
if(filter_mono_exonic) {
  post <- post %>% filter(exons > 1 | filter_pass == FALSE)
  message(" * Removing mono-exonic transcripts")
}

# Only keeping transcripts that pass the filter
filtered_post <- filter(post, filter_pass == TRUE)
message(" * Filtering transcripts based on SQANTI annotation")
message(" * Kept ", nrow(filtered_post), " out of ", nrow(pre), " transcripts (", 
        round(nrow(filtered_post)/nrow(pre) * 100, 1), "%)")


fsm_count <- filtered_post %>% 
  filter(structural_category == "full-splice_match") %>% 
  nrow()

# Calculating percentage of FSM among kept isoforms
fsm_percentage <- round((fsm_count / nrow(filtered_post)) * 100, 1)

# Printing the results
message(" * FSM isoforms kept: ", fsm_count, " out of ", nrow(filtered_post), 
        " total kept isoforms (", fsm_percentage, "%)")


category_counts <- filtered_post %>%
  group_by(structural_category) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Printing the category breakdown
message(" * Structural category breakdown of kept isoforms:")
for(i in 1:nrow(category_counts)) {
  message("   - ", category_counts$structural_category[i], ": ", 
          category_counts$count[i], " (", 
          round(category_counts$count[i]/nrow(filtered_post)*100, 1), "%)")
}

# GENERATING FILTERING REPORT 
# Creating a detailed report of why transcripts were filtered out
reasons <- pre %>%
  filter(!(isoform %in% filtered_post$isoform)) %>%
  transmute(
    isoform, 
    structural_category,
    reasons = case_when(
      structural_category != "full-splice_match" & RTS_stage == TRUE ~ "non-FSM_RT_switching",
      structural_category != "full-splice_match" & diff_to_gene_TTS != 0 & 
        (N_A_downstream_TTS >= 6 | perc_A_downstream_TTS >= 60 | rc_score > 15) ~ "non-FSM_unreliable_TTS",
      exons == 1 & filter_mono_exonic ~ "mono-exonic",
      TRUE ~ "other"
    )
  )

# Summary statistics
message(" * Rejection reasons summary:")
message("   - Non-FSM with RT switching: ", sum(reasons$reasons == "non-FSM_RT_switching"))
message("   - Non-FSM with unreliable TTS: ", sum(reasons$reasons == "non-FSM_unreliable_TTS"))
message("   - Mono-exonic transcripts: ", sum(reasons$reasons == "mono-exonic"))
message("   - Other reasons: ", sum(reasons$reasons == "other"))

mono_exonic_analysis <- pre %>%
  filter(exons == 1) %>%
  group_by(structural_category) %>%
  summarise(
    count = n(),
    percentage = round(n() / sum(pre$exons == 1) * 100, 1)
  ) %>%
  arrange(desc(count))

# Printing % mono-exonic transcirpts
message("Total mono-exonic transcripts: ", sum(pre$exons == 1), 
        " (", round(sum(pre$exons == 1)/nrow(pre)*100, 1), "% of all transcripts)")
message("\nBreakdown by structural category:")
for(i in 1:nrow(mono_exonic_analysis)) {
  message("   - ", mono_exonic_analysis$structural_category[i], ": ", 
          mono_exonic_analysis$count[i], " (", 
          mono_exonic_analysis$percentage[i], "% of mono-exonic)")
}
