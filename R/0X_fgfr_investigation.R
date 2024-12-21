print_banner("FGFR Investigation...")

cna_fgfr2 <- process_cna_gene(data_cna, "FGFR2")
cna_fgfr3 <- process_cna_gene(data_cna, "FGFR3")
cna_fgfr4 <- process_cna_gene(data_cna, "FGFR4")

# Identify ERBB2+ and ERBB2- patients
erbb2_positive <- colnames(cna_erbb2)[cna_erbb2[1, ] > amplification_threshold]
erbb2_negative <- colnames(cna_erbb2)[cna_erbb2[1, ] <= amplification_threshold]

# Check FGFR amplifications in ERBB2+ and ERBB2- patients
fgfr2_amplified_in_erbb2_pos <- cna_fgfr2[1, erbb2_positive] > amplification_threshold
fgfr3_amplified_in_erbb2_pos <- cna_fgfr3[1, erbb2_positive] > amplification_threshold
fgfr4_amplified_in_erbb2_pos <- cna_fgfr4[1, erbb2_positive] > amplification_threshold

fgfr2_amplified_in_erbb2_neg <- cna_fgfr2[1, erbb2_negative] > amplification_threshold
fgfr3_amplified_in_erbb2_neg <- cna_fgfr3[1, erbb2_negative] > amplification_threshold
fgfr4_amplified_in_erbb2_neg <- cna_fgfr4[1, erbb2_negative] > amplification_threshold

# Create contingency tables for FGFR2
fgfr2_table <- table(
  ERBB2 = c(rep("Positive", length(erbb2_positive)), rep("Negative", length(erbb2_negative))),
  FGFR2 = c(fgfr2_amplified_in_erbb2_pos, fgfr2_amplified_in_erbb2_neg)
)
print(fgfr2_table)

# Repeat for FGFR3
fgfr3_table <- table(
  ERBB2 = c(rep("Positive", length(erbb2_positive)), rep("Negative", length(erbb2_negative))),
  FGFR3 = c(fgfr3_amplified_in_erbb2_pos, fgfr3_amplified_in_erbb2_neg)
)
print(fgfr3_table)

# Repeat for FGFR4
fgfr4_table <- table(
  ERBB2 = c(rep("Positive", length(erbb2_positive)), rep("Negative", length(erbb2_negative))),
  FGFR4 = c(fgfr4_amplified_in_erbb2_pos, fgfr4_amplified_in_erbb2_neg)
)
print(fgfr4_table)

# Perform Fisher's Exact Test for FGFR2
fisher_test_fgfr2 <- fisher.test(fgfr2_table)
print(fisher_test_fgfr2)

fisher_test_fgfr3 <- fisher.test(fgfr3_table)
print(fisher_test_fgfr3)

fisher_test_fgfr4 <- fisher.test(fgfr4_table)
print(fisher_test_fgfr4)