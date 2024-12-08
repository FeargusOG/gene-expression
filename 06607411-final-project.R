# Read the RNA-seq file
data_Rnaseq <- read.delim("brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")

# Read the patient data. Skip first 4 rows as they don't contain observations.
data_patient <- read.delim("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4)

# Read the Copy Number Aberrations Data
data_cna <- read.delim("brca_tcga_pan_can_atlas_2018/data_cna.txt")

print(head(data_Rnaseq))
print(head(data_patient))
print(head(data_cna))
