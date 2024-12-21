# 4. Read the Patient Data file: data_clinical_patient.txt
## Skip first 4 rows as they don't contain observations.
data_patient <- read.delim("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4)

# 5. Read the Copy Number Aberrations Data: data_cna.txt
data_cna <- read.delim("brca_tcga_pan_can_atlas_2018/data_cna.txt")
erbb2_cna <- data_cna[data_cna$Hugo_Symbol == "ERBB2", ] |>
  dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id)
rownames(erbb2_cna) <- "ERBB2"

# Transpose the matrix to have patients as rows and genes as columns
transpose_erbb2_cna <- t(erbb2_cna)

# Standardise IDs in the transpose_erbb2_cna dataset
transpose_erbb2_cna_rownames <- rownames(transpose_erbb2_cna) # Extract row names (patient IDs)
transpose_erbb2_cna_rownames <- sub("\\.\\d+$", "", transpose_erbb2_cna_rownames) # Remove trailing ".01"
transpose_erbb2_cna_rownames <- gsub("\\.", "-", transpose_erbb2_cna_rownames) # Replace '.' with '-'

# Update rownames in the transpose_erbb2_cna object
rownames(transpose_erbb2_cna) <- transpose_erbb2_cna_rownames


# Set the PATIENT_ID column as the rownames
Y_patient <- data.frame(
  patient_id = data_patient$PATIENT_ID,
  time = data_patient$OS_MONTHS,
  status = data_patient$OS_STATUS
)
rownames(Y_patient) <- Y_patient$patient_id

# Drop the patient_id column
Y_patient$patient_id <- NULL
# Ensure the OS_MONTHS column is numeric and > 0
Y_patient$time <- as.numeric(Y_patient$time)
Y_patient <- Y_patient[Y_patient$time > 0,]
# Conver the OS_STATUS to numeric 1 or 0
Y_patient$status <- ifelse(Y_patient$status == "1:DECEASED", 1, 
                           ifelse(Y_patient$status == "0:LIVING", 0, NA))



# Merge Y_patient and transpose_erbb2_cna by patient ID (rownames)
Y_erbb2_amp <- merge(Y_patient, transpose_erbb2_cna, by = "row.names", all = FALSE)

# Set the merged rownames to match the patient IDs
rownames(Y_erbb2_amp) <- Y_erbb2_amp$Row.names
Y_erbb2_amp$Row.names <- NULL

# Ensure time and status columns are numeric
Y_erbb2_amp$time <- as.numeric(Y_erbb2_amp$time)  # Ensure time is numeric
Y_erbb2_amp$status <- as.numeric(Y_erbb2_amp$status)  # Ensure status is numeric (1 = event, 0 = censored)

# Create the erbb2_amp column as a binary categorical variable
dat_test_erbb2 <- data.frame(
  time = Y_erbb2_amp[, "time"], 
  status = Y_erbb2_amp[, "status"], 
  erbb2_amp = ifelse(Y_erbb2_amp[, "ERBB2"] > 0, "ERBB2+", "ERBB2-")  # Create binary column
)

# Fit a survival curve using the test data
s_fit <- survfit(Surv(time, status) ~ erbb2_amp, data = dat_test_erbb2)

# Plot the survival curve
ggsurvplot(
  s_fit, 
  data = dat_test_erbb2, 
  # conf.int = TRUE, 
  # risk.table = TRUE, 
  title = "Kaplan-Meier Survival by ERBB2 Amplification",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  legend.labs = c("ERBB2-", "ERBB2+"), # Customise legend labels
  palette = c("darkblue", "red")       # Customise colours
)