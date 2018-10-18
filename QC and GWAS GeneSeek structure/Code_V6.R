#----------- Libraries ---------#
library (GWASTools)
library(ade4)
library(SNPRelate)
#-------------------------------#


#--------- Data paths ----------#
# Path where are stored the genotypes. 
path.folders = "../Data"

# File with genetic origin and phenotypes
Scan_Annotation_Data_Frame = read.csv(file= "../Data/Genetic origin and continuous phenotype - Toy Data.csv", header = TRUE, sep=",")

# Folder for outputs
path_outputs <- "../Results"
dir.create(path_outputs)
#-------------------------------#


#--------- Parameters ----------#
# The following parameters are mandatory

# Number of autosomes of the specie
autosomes <- 29

# Total number of SNPs in autosomes
num_SNP_auto = 29794

# Quality control
## Minor allele frequency curoff
maf = 0.01
## Hardy-Weinberg equilibrium test. P-value cutoff 
pvalue = 0.05
## GenCall Score cutooff
mean_GC.score_sample = 0.8 # GenCall Score mean per sample
median_GC.score_sample = 0.865 # GenCall Score median per sample
mean_GC.score_snp = 0.8 # GenCall Score mean per SNP
median_GC.score_snp = 0.7 # GenCall Score median per SNP
## Cutoff of missing call rate
cutoof=0.05
## Significance of population stratification's manhattan plot
pop_signif = 1e-7

# GWAS parameters
## Phenotype column name
outcome = "trait"
## Model type.
model.type = "linear"
## Covariates.
covar = NULL
## Covariate with interaction with genotype. It must be a name of covar
ivar = NULL
## Confidence interval
CI = 0.95
## Number of SNPs to read in at once
block.size = 5000
## Multiple comparision test. fdr or bonferroni
method = "fdr"
## Manhattan plot significance
signif = 1e-5 
#-----------------------------#


#----------- Uncompress genotype files -----------#

# Create the list of batch folders
Main_folders.zip = list.files(path.folders, pattern="/*.zip", recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
# Uncompress batch folders
lapply(Main_folders.zip, 
       function(Main_folders) unzip ( Main_folders , exdir = path.folders))
# List of attached files
files.zip = list.files(path.folders, pattern="/*.zip", recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
# Extract files
lapply(files.zip, 
       function(files) unzip ( files , exdir = path.folders))
#------------------------------------------------#


#----------- Genotype files input ---------------#

# LocusSummary
## Function
read_function_LocusSummary = function (file, rows) {
  if (missing(rows)) {
    read.csv(file=file, skip=2, header = TRUE, sep=",")} 
  else {
    read.csv(file=file, skip=rows, header = TRUE, sep=",")}}
## Files list
File_names_LocusSummary = list.files(path.folders, pattern="/*LocusSummary.csv", recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
## Input files
LocusSummary = lapply(File_names_LocusSummary, function(file) read_function_LocusSummary(file))

# SNP Map
## Function
read_function_SNP_Map = function (file, rows) { 
  if (missing(rows)) {
    read.table(file=file, sep = "\t" , header = TRUE, check.names = FALSE, blank.lines.skip = TRUE, na.strings = "NA")} 
  else {
    read.table(file=file, sep = "\t" , header = TRUE, check.names = FALSE, blank.lines.skip = TRUE, na.strings = "NA", skip=rows)}}
## Files list
File_names_SNP_Map = list.files(path.folders, pattern="/*SNP_Map.txt", recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
## Input files
SNP_Map = lapply(File_names_SNP_Map, function(file) read_function_SNP_Map(file)) 

# FinalReport
## Function
read_function_FinalReport = function (file, rows) { 
  if (missing(rows)) {
    read.table(file=file, sep = "\t" , header = TRUE, check.names = FALSE, blank.lines.skip = TRUE, na.strings = "NA", skip=9)} 
  else {
    read.table(file=file, sep = "\t" , header = TRUE, check.names = FALSE, blank.lines.skip = TRUE, na.strings = "NA", skip=rows)}}
## Files list
File_names_FinalReport = list.files(path.folders, pattern="/*FinalReport.txt", recursive = TRUE, full.names = TRUE, include.dirs = FALSE) 
## Input files
FinalReport <- lapply(File_names_FinalReport, function(file) read_function_FinalReport(file))
## Merge FinalReports files
FinalReport <- do.call("rbind", FinalReport)
## Remove SNP that are not in all the samples
### frequency of SNP Name
freq_SNP <- as.data.frame(table(FinalReport$"SNP Name")) 
### list of SNPs in all frequency
SNPs_in_all <- freq_SNP[freq_SNP$Freq >= length(unique(FinalReport$"Sample ID")),]
### Filter FinalReport
FinalReport <- FinalReport[FinalReport$"SNP Name" %in% SNPs_in_all$Var1,]
### Export and import the information to avoid error...
write.table(FinalReport, file = "FinalReport.txt")
FinalReport <- read.table("FinalReport.txt", header = TRUE, na.strings = TRUE)

# FinalReportCNV
## Function
read_function_FinalReportCNV = function (file, rows) { 
  if (missing(rows)) {
    read.csv(file=file, skip=9, header = TRUE, sep=",")} 
  else {
    read.csv(file=file, skip=rows, header = TRUE, sep=",")}}
## Files list
File_names_FinalReportCNV = list.files(path.folders, pattern="/*FinalReportCNV.csv", recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
## Input files
FinalReportCNV = lapply(File_names_FinalReportCNV, function(file) read_function_FinalReportCNV(file))
## Merge all FinalReportsCNV
FinalReportCNV <- do.call("rbind", FinalReportCNV)
## Remove SNP that are not in all the samples
### Filter FinalReportCNV
FinalReportCNV <- FinalReportCNV[FinalReportCNV$"SNP.Name" %in% SNPs_in_all$Var1,]
### Export and import the information to avoid error...
write.table(FinalReportCNV, file = "FinalReportCNV.txt")
FinalReportCNV <- read.table("FinalReportCNV.txt", header = TRUE, na.strings = TRUE)
#------------------------------------------------#


#-------- Creation of GWASTools objects ---------#

# SNP Annotation Data
merge_SNP_Map_and_LocusSummary <- merge(x = SNP_Map[[1]], y = LocusSummary[[1]], by.x = "Name", by.y = "Locus_Name")
## Rename chromosomes
Chromosome.vector <- merge_SNP_Map_and_LocusSummary[3]
Chromosome.vector <- Chromosome.vector [,]
Chromosome.vector <- as.character(Chromosome.vector)
Chromosome.vector <- replace(Chromosome.vector,which(Chromosome.vector=="X"),autosomes+1)
Chromosome.vector <- replace(Chromosome.vector,which(Chromosome.vector=="XY"),autosomes+2)
Chromosome.vector <- replace(Chromosome.vector,which(Chromosome.vector=="Y"),autosomes+3)
Chromosome.vector <- replace(Chromosome.vector,which(Chromosome.vector=="MT"),autosomes+4)
Chromosome.vector <- replace(Chromosome.vector,which(!Chromosome.vector %in% 1:autosomes+4),autosomes+5)
## Creation of the snp anottation data frame
snp_annot_data_frame <- data.frame(
  snpName = merge_SNP_Map_and_LocusSummary$Name,
  chromosome = as.integer (Chromosome.vector),
  position = merge_SNP_Map_and_LocusSummary$Position,
  alleleA = as.character(substr(merge_SNP_Map_and_LocusSummary$SNP, 2, 2)),
  alleleB = as.character(substr(merge_SNP_Map_and_LocusSummary$SNP, 4, 4)),
  BeadSetID = as.integer(merge_SNP_Map_and_LocusSummary$Illumicode_Name),
  tAA = merge_SNP_Map_and_LocusSummary$AA_T_Mean,
  tAB = merge_SNP_Map_and_LocusSummary$AB_T_Mean,
  tBB = merge_SNP_Map_and_LocusSummary$BB_T_Mean,
  rAA = merge_SNP_Map_and_LocusSummary$AA_R_Mean,
  rAB = merge_SNP_Map_and_LocusSummary$AB_R_Mean,
  rBB = merge_SNP_Map_and_LocusSummary$BB_R_Mean)
## Remove SNP that are not in all the samples
### Filter FinalReport
snp_annot_data_frame <- snp_annot_data_frame[snp_annot_data_frame$"snpName" %in% SNPs_in_all$Var1,]
### Export and import the information to avoid error...
write.table(snp_annot_data_frame, file = "snp_annot_data_frame.txt")
snp_annot_data_frame <- read.table("snp_annot_data_frame.txt", header = TRUE, na.strings = TRUE)
## Creation of snpID as an integer sorted by chromosome and position
### Sort snp_annot_data_frame by chromosome and then by position
snp_annot_data_frame <- snp_annot_data_frame[order( snp_annot_data_frame$chromosome, snp_annot_data_frame$position ),]
### Cration of snpID
snpID = as.integer(c(1:nrow(snp_annot_data_frame)))
## Merge snpID whit snp_annot_data_frame
snp_annot_data_frame_snpID <- data.frame(
  snpID = snpID,
  snp_annot_data_frame)
##  Creation of the SNP Annotation Data Object.
autosomes <- as.integer(autosomes)
snpAnnot <- SnpAnnotationDataFrame(
  snp_annot_data_frame_snpID, 
  autosomeCode = 1:autosomes,
  XchromCode = as.integer(autosomes+1),
  XYchromCode = as.integer(autosomes+2),
  YchromCode = as.integer(autosomes+3),
  MchromCode = as.integer(autosomes+4))
## Add metadata to describe the columns snpAnnot.
meta <- varMetadata(snpAnnot)
meta[c("snpID", "snpName", "chromosome", "position", "alleleA", "alleleB",
       "BeadSetID", "tAA", "tAB", "tBB", "rAA", "rAB", "rBB"),
     "labelDescription"] <- c("unique integer ID for SNPs, is the Index column of SNP_Map",
                              "snpName of Illumina",
                              paste("integer code for chromosome: 1:29=autosomes,",
                                    "30=X, 31=pseudoautosomal, 32=Y, 33=Mitochondrial, 34=Unknown"), 
                              "base pair position on chromosome",
                              "alelleA", "alleleB",
                              "BeadSet ID from Illumina",
                              "mean theta for AA cluster",
                              "mean theta for AB cluster",
                              "mean theta for BB cluster",
                              "mean R for AA cluster",
                              "mean R for AB cluster",
                              "mean R for BB cluster")
varMetadata(snpAnnot) <- meta

# Scan Annotation Data
## Creation of scanID
scanID <- as.numeric(c(row.names(Scan_Annotation_Data_Frame)))
## Merge scanID with the Scan_Annotation_Data_Frame
Scan_Annotation_Data_Frame_scanID <- data.frame (scanID , Scan_Annotation_Data_Frame)
## Export and import scan.annotation to avoid error...
write.table(Scan_Annotation_Data_Frame_scanID, file = "Scan_Annotation_Data_Frame_scanID.txt") 
Scan_Annotation_Data_Frame_scanID  <- read.table("Scan_Annotation_Data_Frame_scanID.txt", header = TRUE, na.strings = TRUE)
## Create the ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(Scan_Annotation_Data_Frame_scanID)
## Add metadata to describe the columsns
meta <- varMetadata(scanAnnot)
meta[c("scanID", "subjectID", "genetic_origin", "sex", "batch", "trait"),
     "labelDescription"] <-
  c("unique ID for scans",
    "subject identifier (may have multiple scans)",
    "genetic origin group",
    "sex coded as M=male and F=female",
    "genotyping batch",
    "trait")
varMetadata(scanAnnot) <- meta

## Creation of Data Files
### Creation of a file for each samples
#### Define the path
path_GenotypeData <- paste(path.folders,"/GenotypeData", sep="") 
dir.create(path_GenotypeData)
#### Data frame of scanID whit Sample.ID
scanID_Sample.ID <- Scan_Annotation_Data_Frame_scanID [1:2]
#### Data frame for each sample of FinalReport
FinalReport_individual<-split(FinalReport, with(FinalReport, FinalReport[2])) 
#### Data frame for each sample of FinalReportCNV
FinalReportCNV_individual<-split(FinalReportCNV, with(FinalReportCNV, Sample.ID))
#### Merge FinalReport whit FinalReportCNV 
func_Merg_FinalReport_FinalReportCNV <- function(x,y){merge(x, y, by = c("Sample.ID", "SNP.Name"))}
FinalReports <- mapply(func_Merg_FinalReport_FinalReportCNV, FinalReport_individual, FinalReportCNV_individual, SIMPLIFY=FALSE)
#### Merge FinalReports whit scanID
func_Merg_FinalReport_scanID <- function(x,y){merge(x, y, by.x="Sample.ID", by.y="subjectID")}
FinalReport_scanID <- lapply(FinalReports, func_Merg_FinalReport_scanID, scanID_Sample.ID)
#### Create files names
func_create_file_names <- function (x) {paste(path_GenotypeData,"/", x[1,"scanID"],".csv", sep="")}
raw_data_file <- lapply(FinalReport_scanID, func_create_file_names)
#### Creation of files
func_creat_files <- function (x,y) {write.csv(x, file = y)}
mapply (func_creat_files, FinalReport_scanID, raw_data_file, SIMPLIFY=FALSE)
#### Inputing and exporting files to avoid (error) quotes in numbers when "scan" is used to create "diag.qxy"
##### Function
func_individual_files = function (file) {
  read.csv(file=file, colClasses = c("factor", "factor", "factor", NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "factor"))}
##### List of file paths
individual_files_paths <- list.files(path = paste(path.folders, "/GenotypeData", sep=""), pattern="*.csv", recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
##### Creation of files
individual_files <- lapply(individual_files_paths, func_individual_files) 
##### File names
individual_files_names <- lapply(individual_files, func_create_file_names)
##### Creation of files
mapply (func_creat_files, individual_files, individual_files_names, SIMPLIFY=FALSE)
### Create scan.annotation
scan.annotation	 <- data.frame(scanID = scanID_Sample.ID$scanID,
                               scanName = scanID_Sample.ID$subjectID,
                               file = paste(scanID_Sample.ID$scanID,".csv", sep=""))
scan.annotation <- scan.annotation[order(scan.annotation$scanID),] 
### Create snp.annotation
snp.annotation <- snp_annot_data_frame_snpID[,c("snpID", "snpName", "chromosome", "position")]
snp.annotation <- snp.annotation[order(snp.annotation$snpID),] 

### Genotype File
geno.file <- "tmp.geno.gds"
col.nums <- as.integer(c(3,4,9,10))
names(col.nums) <- c("sample", "snp", "a1", "a2")
diag.geno.file <- "diag.geno.RData" # Name of the output file
diag.geno <- createDataFile (path = path_GenotypeData,
                             filename = geno.file,
                             file.type= "gds",
                             variables = "genotype",
                             snp.annotation = snp.annotation, 
                             scan.annotation = scan.annotation,
                             sep.type = ",",
                             skip.num = 1,
                             col.total = 16,
                             col.nums = col.nums,
                             scan.name.in.file = 0,
                             diagnostics.filename = diag.geno.file,
                             verbose = FALSE)

### Intensity File
qxy.file <- "tmp.qxy.gds"
col.nums <- as.integer(c(3,4,11,12,13))
names(col.nums) <- c("sample", "snp", "quality", "X", "Y")
diag.qxy.file <- "diag.qxy.RData"
diag.qxy <- createDataFile(path = path_GenotypeData,
                           filename = qxy.file,
                           file.type = "gds",
                           variables = c("quality","X","Y"),
                           snp.annotation = snp.annotation,
                           scan.annotation = scan.annotation,
                           sep.type = ",",
                           skip.num = 1,
                           col.total = 16, 
                           col.nums = col.nums,
                           scan.name.in.file = 0,
                           diagnostics.filename = diag.qxy.file,
                           verbose = FALSE)

### BAF and LRR file
bl.file <- "tmp.bl.gds"
col.nums <- as.integer(c(3,4,14,15))
names(col.nums) <- c("sample", "snp", "BAlleleFreq", "LogRRatio")
diag.bl.file <- "diag.bl.RData"
diag.bl <- createDataFile(path = path_GenotypeData,
                          filename = bl.file,
                          file.type = "gds",
                          variables = c("BAlleleFreq","LogRRatio"),
                          snp.annotation = snp.annotation,
                          scan.annotation = scan.annotation,
                          sep.type = ",",
                          skip.num = 1,
                          col.total = 16,
                          col.nums = col.nums,
                          scan.name.in.file = 0,
                          diagnostics.filename = diag.bl.file,
                          verbose = FALSE)

### Creation of genoData.
#### Initialize GdsGenotypeReader file
(gds <- GdsGenotypeReader(geno.file))
#### Renumber chromosomes
gds@autosomeCode <- as.integer(c(1:autosomes))
gds@XchromCode <- as.integer(autosomes+1)
gds@XYchromCode <- as.integer(autosomes+2)
gds@YchromCode <- as.integer(autosomes+3)
gds@MchromCode <- as.integer(autosomes+4)
#### Creat genoData
genoData <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
### Create first IntensityData. It is a combination of IntensityReader with SNP and Scan annotation
#### GdsIntensityReader file
(gds <- GdsIntensityReader(qxy.file))
#### Renumber chromosomes
gds@autosomeCode <- as.integer(c(1:autosomes))
gds@XchromCode <- as.integer(autosomes+1)
gds@XYchromCode <- as.integer(autosomes+2)
gds@YchromCode <- as.integer(autosomes+3)
gds@MchromCode <- as.integer(autosomes+4)
#### Create first IntensityData
qxyData <- IntensityData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
### Create second Intensity file
#### GdsIntensityReader file
(gds <- GdsIntensityReader(bl.file))
#### Renumber chromosomes
gds@autosomeCode <- as.integer(c(1:autosomes))
gds@XchromCode <- as.integer(autosomes+1)
gds@XYchromCode <- as.integer(autosomes+2)
gds@YchromCode <- as.integer(autosomes+3)
gds@MchromCode <- as.integer(autosomes+4)
blData <- IntensityData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
#---------------------------------------


#---------- Quality control ------------
# Quality control path
path_quality_control <- paste(path_outputs,"/Quality control", sep="") 
dir.create(path_quality_control)


# SNP quality control

## Create the results path
path_SNP_quality <- paste(path_quality_control,"/SNP quality control", sep="") 
dir.create(path_SNP_quality)

## GC score analysis per SNP
### Create the results path
SNP_genotype_quality_scores <- paste(path_SNP_quality,"/SNP GenCall Score", sep="") 
dir.create(SNP_genotype_quality_scores)
### Calculate mean and median GC score per SNP
qual.results_snp <- data.frame(qualityScoreBySnp(qxyData, genoData))
### Create sample column
qual.results_snp <- data.frame(snpID = rownames(qual.results_snp),
                               qual.results_snp)
### Filter out SNPs with mean and median bellow cutoff
qual.results_snp_filt <- qual.results_snp[qual.results_snp$mean.quality >= mean_GC.score_snp,]
qual.results_snp_filt <- qual.results_snp_filt[
  qual.results_snp_filt$median.quality >= median_GC.score_snp,]
### Export list of SNPs with mean and median greater than cutoff
write.table(format(qual.results_snp_filt, digits = 2), 
            file = paste (SNP_genotype_quality_scores,
                          "/1. SNPs with mean and median greater than cutoff.txt", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE) 
### Export plots of number of samples before and after the filter
set.seed(1)
pdf(paste(SNP_genotype_quality_scores, 
          "/02. Number of SNPs before and after the filter.pdf", sep=""), 7, 5)
par(mfrow=c(2,2))
hist(qual.results_snp$mean.quality, 
     main="Mean GC score of SNP before filter",
     xlab="Mean GC score",
     ylab = "Number of SNPs")
hist(qual.results_snp_filt$mean.quality, 
     main="Mean GC score of SNPs after filter",
     xlab="Mean GC score",
     ylab = "Number of SNPs")
hist(qual.results_snp$median.quality, 
     main="Median GC score of SNPs before filter",
     xlab="Median GC score",
     ylab = "Number of SNPs")
hist(qual.results_snp_filt$median.quality, 
     main="Median GC score of SNPs after filter",
     xlab="Median GC score",
     ylab = "Number of SNPs")
dev.off()


## Population stratification analysis
### Create folders for results
path_stratification_analysis <- paste(path_SNP_quality,"/Population stratification", sep="") 
dir.create(path_stratification_analysis)
### Multiple correspondence analysis before filter
#### Create de data frame
data_concat <- data.frame(SNP.Name = FinalReport$SNP.Name,
                          Sample.ID = FinalReport$Sample.ID,
                          Allele1...AB = FinalReport$Allele1...AB,
                          Allele2...AB = FinalReport$Allele2...AB)
#### Unify columns of genotypes
data_concat <- data.frame(data_concat[,1:2], 
                          Genotype = paste(data_concat[,3], data_concat[,4], sep=''))
#### Transform genotype from AA, AB and BB format to  0, 1 and 2, respectively.
data.ch <- as.character(data_concat[, 3])
Genotype <- replace(data.ch,which(data.ch=="AA"),2)
Genotype <- replace(Genotype,which(data.ch=="AB"),1)
Genotype <- replace(Genotype,which(data.ch=="BB"),0)
Genotype <- replace(Genotype,which(data.ch=="--"),NA)
#### Create a data frame with new genotype format
data_210 <- data.frame(data_concat[1:2], Genotype = Genotype)
#### Transpose database: SNPs in columns and individuals in rows
data_split <- data.frame(subjectID = unique(data_210$Sample.ID), 
                         split(data_210$Genotype, f = data_210$SNP.Name))
#### Export and import
write.table(data_split, file = "data_split.txt") 
data_split <- read.table("data_split.txt", header = TRUE, na.strings = TRUE)
#### Merge genotype whit genetic_origin
data_origin <- merge (x = Scan_Annotation_Data_Frame_scanID[2:3], y = data_split, by = "subjectID") 
#### Remove unnecessary columns
data_origin <- data_origin [,-c(1)] 
#### Export and import
write.table(data_origin, file = "data_origin.txt")
data_origin_imp <- read.table("data_origin.txt", header = TRUE, na.strings = TRUE)
#### Transform in numbers in factors
data_origin_imp2 <- data.frame(apply(data_origin_imp, 2, as.factor))
#### Make the multiple correspondence analysis
acm_origin <- dudi.acm(data_origin_imp2, row.w = rep(1, nrow(data_origin_imp2)), scan = FALSE, nf =2)
#### Generate the scatter plot of before filter data
set.seed(1)
pdf(paste (path_stratification_analysis,
           "/01. MCA_before_filter.pdf", sep=""), 7, 5)
s.class(acm_origin$li,data_origin_imp2$genetic_origin)
dev.off()
### Remove SNPs associated to genetic origin
#### Association analysis
Aso_Orig<-assocRegression(genoData,
                          outcome = "genetic_origin",
                          model.type="logistic",
                          snpStart = 1,
                          snpEnd = num_SNP_auto)
#### Manhattan Plot before filter
set.seed(1)
pdf(paste (path_stratification_analysis, "/02. Manhattan_Plot_beforer_filter.pdf", sep=""), 7, 5)
manhattanPlot(Aso_Orig$Wald.pval, Aso_Orig$chr, signif=pop_signif)
dev.off()
#### Remove SNPs associated to genetic origin
Aso_Orig_sort <- sort(Aso_Orig$Wald.pval)
SNP_origin <- Aso_Orig [Aso_Orig$Wald.pval <= pop_signif,] 
SNP_origin <- SNP_origin[!is.na(SNP_origin)]
SNP_no_origin <- Aso_Orig [Aso_Orig$Wald.pval > pop_signif,]
#### Manhattan Plot after filter
set.seed(1)
pdf(paste (path_stratification_analysis, "/03. Manhattan_Plot_after_filter.pdf", sep=""), 7, 5)
manhattanPlot(SNP_no_origin$Wald.pval, SNP_no_origin$chr, signif=pop_signif)
dev.off()
### Scatter plot after remove SNPs associated to origin
#### Remove columns of snpID. 		
SNPs_to_remove <- unique(SNP_origin$snpID) 
SNPs_to_remove <- data.frame ("snpID" = SNPs_to_remove [ !is.na(SNPs_to_remove)] )
SNPs_to_remove <- merge (x = SNPs_to_remove, y = snp_annot_data_frame_snpID [1:2], by = "snpID")
#### Export list of removed SNPs
write.table(SNPs_to_remove, 
            file = paste (path_stratification_analysis, "/04. SNP_associated_with_origin.txt", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)
#### Identify columns with SNPs to remove
#### Transform names to match with columns of "data_origin_imp2"
name_SNPs_to_remove <- gsub ("-", ".", SNPs_to_remove[,2])
#### Identify the number of columns of SNPs to remove
col_SNPs_to_remove <-  match(name_SNPs_to_remove, names(data_origin_imp2))
#### Remove columns with SNPs associated to genetic origin
data_snp_origin <- data_origin_imp2 [ , -col_SNPs_to_remove]
#### Stratification analysis after remove SNPs
acm_origin <- dudi.acm(data_snp_origin, row.w = rep(1, nrow(data_snp_origin)), scan = FALSE, nf =2)
#### print and store
set.seed(1)
pdf(paste (path_stratification_analysis, "/05. Stratification_analysis_after_filter.pdf", sep=""), 7, 5)
s.class(acm_origin$li,data_origin_imp2$genetic_origin)
dev.off()

## Hardy-Weinberg equilibrium
### Create result folder
HW_equilibrium <- paste(path_SNP_quality,"/HW_equilibrium", sep="") 
dir.create(HW_equilibrium)
###  Separate founders individuals
nonfounders <- scanAnnot$sire != 0 & scanAnnot$dam != 0
### Chromosome vector
chr <- getChromosome(genoData)
### Autosome vector
auto <- range(which(chr %in% 1:autosomes))
### X chromosome vector
X <- range(which(chr == autosomes+1))
### Hardy-Weinberg test for autosomes
hwe <- exactHWE(genoData, snpStart=auto[1], snpEnd=auto[2])
### Hardy-Weinberg test for X chromosome
hweX <- exactHWE(genoData, snpStart=X[1], snpEnd=X[2])
### Merge tests of autosomes and X chromosome
hwe <- rbind(hwe, hweX)
### Create column with the total of individuals per SNP analyzed
hwe$N <- hwe$nAA + hwe$nAB + hwe$nBB
### Change column names
colnames(hwe) <- c("snpID", "chromosome", "individuals_AA", "individuals_AB", "individuals_BB", "Minor_allele_frequency",
                   "menor.allele", "Inbreeding_coefficient", "p-val", "total_individuals")
### Export Hardy-Weinberg before filter equilibrium
write.table(hwe, file = paste (HW_equilibrium, "/01. Hardy-Weinberg equilibrium before filter.txt", sep=""),
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE)
### Filter data with low p-value
hwe_filtered <- hwe[hwe$`p-val` >= pvalue,]
### Export histogram of inbreeding coefficient
set.seed(1)
pdf(paste (HW_equilibrium, "/02. Inbreeding coefficient.pdf", sep=""), 7, 5)
hist(hwe$Inbreeding_coefficient, main="Inbreeding coefficient before filter",
     xlab="Inbreeding Coefficient")
hist(hwe_filtered$Inbreeding_coefficient, main="Inbreeding coefficient after filter",
     xlab="Inbreeding Coefficient")
dev.off()
### Export Hardy-Weinberg equilibrium after filter
write.table(hwe_filtered, file = paste (HW_equilibrium, "/03. Hardy-Weinberg equilibrium after filter.txt", sep=""),
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

## Menor allele frequency analysis
### Create results path
path_MAF_analysis <- paste(path_SNP_quality,"/Menor allele frequency analysis", sep="") 
dir.create(path_MAF_analysis)
### Expor MAF before filter
set.seed(1)
pdf(paste(path_MAF_analysis, "/01. Minor allele frequency before filter.pdf", sep=""), 7, 5)
hist(hwe$Minor_allele_frequency, nclass = 20, main = "Minor allele frequency before filter",
     xlab = "Minor allele frequency", ylab = "Number of SNPs" )
dev.off()
### Filter out SNPs with MAF less than the cutoof
hwe_filtered_MAF <- hwe[hwe$Minor_allele_frequency>=maf,]
#### Export SNPs with MAF greater than cutoff 
write.table(hwe[hwe$Minor_allele_frequency<maf,-c(8:9)], 
            file = paste (path_MAF_analysis,
                          "/02. SNPs with MAF greater than cutoff", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE) 
### Expor MAF after filter
set.seed(1)
pdf(paste(path_MAF_analysis, "/03. Minor allele frequency after filter.pdf", sep=""), 7, 5)
hist(hwe_filtered_MAF$Minor_allele_frequency,  nclass = 20, 
     main = "Minor allele frequency after filter",
     xlab = "Minor allele frequency", ylab = "Number of SNPs" )
dev.off()


# Sample Quality Checks
## Folders to store the outpusts
Sample_Quality_Checks <- paste(path_quality_control,"/Sample_Quality_Checks", sep="") 
dir.create(Sample_Quality_Checks)
Sample_genotype_quality_scores <- 
  paste(Sample_Quality_Checks,"/Sample_genotype_quality_scores", sep="") 
dir.create(Sample_genotype_quality_scores)
B_Allele_Frequency_variance_analysis <-
  paste(Sample_Quality_Checks,"/B_Allele_Frequency_variance_analysis", sep="") 
dir.create(B_Allele_Frequency_variance_analysis)
sample_missingness_and_heterozygosities <-
  paste(Sample_Quality_Checks,"/sample_missingness_and_heterozygosities", sep="") 
dir.create(sample_missingness_and_heterozygosities)

## GC score per sample
### Calculate mean and median GC score per sample
qual.results <- data.frame(qualityScoreByScan(qxyData, genoData))
### Create sample column
qual.results <- data.frame(Sample = rownames(qual.results),
                           qual.results)
### Filter out samples with mean and median bellow cutoff
qual.results_filt <- qual.results[qual.results$mean.quality >= mean_GC.score_sample,]
qual.results_filt <- qual.results_filt[qual.results_filt$median.quality >= median_GC.score_sample,]
### Export list of samples with mean and median greater than cutoff
write.table(format(qual.results_filt, digits = 2), 
            file = paste (Sample_genotype_quality_scores,
                          "/1. Samples with mean and median greater than cutoff.txt", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE) 
### Export plots of number of samples before and after the filter
set.seed(1)
pdf(paste(Sample_genotype_quality_scores, 
           "/02. Number of samples before and after the filter.pdf", sep=""), 7, 5)
par(mfrow=c(2,2))
hist(qual.results$mean.quality, 
     main="Mean GC score of Samples before filter",
     xlab="Mean GC score",
     ylab = "Number of samples")
hist(qual.results_filt$mean.quality, 
     main="Mean GC score of Samples after filter",
     xlab="Mean GC score",
     ylab = "Number of samples")
hist(qual.results$median.quality, 
     main="Median GC score of Samples before filter",
     xlab="Median GC score",
     ylab = "Number of samples")
hist(qual.results_filt$median.quality, 
     main="Median GC score of Samples after filter",
     xlab="Median GC score",
     ylab = "Number of samples")
dev.off()

## B Allele Frequency variance analysis
### Number of bins for each chromosome. Each chromosome is divided into 12 sections with equal numbers of SNPs.
nbins <- rep(12, autosomes+1) 
### Calculate the standard deviation of BAF at each window in each sample
slidingBAF12 <- sdByScanChromWindow(blData, genoData, nbins=nbins)
### Calculate the mean of BAF at each window in each sample
sds.chr <- meanSdByChromWindow(slidingBAF12, scanAnnot$sex)
### Identify windows with BAF standard deviations four time higher compared to other samples. 
res12bin4sd <- findBAFvariance(sds.chr, slidingBAF12, scanAnnot$sex, sd.threshold=3)
### Export a list unusually high BAF standard deviation sample-chromosome windows
unusually_high_BAF_sd <- data.frame (res12bin4sd)
write.table(unusually_high_BAF_sd, 
            file = paste (B_Allele_Frequency_variance_analysis,
                          "/01. Unusually high BAF standard deviation sample-chromosome", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE) 
### Plot unusually high BAF standard deviation sample-chromosome windows.
scanID <- as.integer(res12bin4sd[, "scanID"])
chrom <- as.integer(res12bin4sd[, "chromosome"])
chrom[res12bin4sd[, "chromosome"] == "X"] <- autosomes+1
bincode <- paste("Bin", res12bin4sd[, "bin"], sep = " ")
set.seed(1)
pdf(paste (B_Allele_Frequency_variance_analysis, 
           "/02. Unusually high BAF standard deviation sample-chromosome.pdf",
           sep=""), 7, 5)
chromIntensityPlot(blData, scanID, chrom, info=bincode, ideogram=FALSE)
dev.off()


## Missingness and heterozygosity within samples. Samples with high heterozygosity may indicate a mixed sample. Also, it can be identified any outliers with regard to missingness.

### Calculate the missing call rate by sample and chromosome.
miss <- missingGenotypeByScanChrom(genoData)
miss.rate <- t(apply(miss$missing.counts, 1, function(x) {x / miss$snps.per.chr}))
miss.rate <- as.data.frame(miss.rate)
cols <- names(miss.rate) %in% c(1:autosomes, "X", "XY", "Y", "M", "U")

### Export the missingness by chromosome
set.seed(1)
pdf(paste (sample_missingness_and_heterozygosities,
           "/01. Missingness by chromosome.pdf", sep=""), 7, 5)
boxplot(miss.rate[,cols], main="Missingness by Chromosome", ylab="Proportion Missing", 
        xlab="Chromosome")
dev.off()

### Export X chromosome missingness by sex
set.seed(1)
pdf(paste (sample_missingness_and_heterozygosities, 
           "/02. X Chromosome Missingness by Sex.pdf", sep=""), 7, 5)
boxplot(miss.rate$X ~ scanAnnot$sex, main="X Chromosome Missingness by Sex", ylab="Proportion Missing")
dev.off()

### Calculate heterozygosity by sample and chromosome
het.results <- hetByScanChrom(genoData)
#### Autosomal and X chromosome heterozygosity. High heterozygosity outliers could be due to sample contamination, while low heterozygosity may due to chromosomal anomalies.
scanAnnot$het.A <- het.results[,"A"]
scanAnnot$het.X <- het.results[,"X"]
#### Add variables description
varMetadata(scanAnnot)["het.A", "labelDescription"] <-
  "fraction of heterozygotes for autosomal SNPs"
varMetadata(scanAnnot)["het.X", "labelDescription"] <-
  "fraction of heterozygotes for X chromosome SNPs"
#### Export autosomal heterozygosity
set.seed(1)
pdf(paste (sample_missingness_and_heterozygosities, "/03. Autosomal Heterozygosity.pdf", sep=""), 7, 5)
boxplot(scanAnnot$het.A ~ scanAnnot$genetic_origin, main="Autosomal Heterozygosity")
dev.off()
#### Export X chromosome heterozygosity in females
set.seed(1)
pdf(paste (sample_missingness_and_heterozygosities, "/04. X Chromosome Heterozygosity in Females.pdf", sep=""), 7, 5)
female <- scanAnnot$sex == "F"
boxplot(scanAnnot$het.X[female] ~ scanAnnot$genetic_origin[female], main="X Chromosome Heterozygosity in Females")
dev.off()


# Sample Identity Checks

## Create result folders
Sample_Identity_Checks <- paste(path_quality_control,"/Sample_Identity_Checks", sep="") 
dir.create(Sample_Identity_Checks)
Misannotated_Sex_Check <- paste(Sample_Identity_Checks,"/Mis-annotated_Sex_Check", sep="") 
dir.create(Misannotated_Sex_Check)

## Determine discrepancies between the annotated sex and genetic sex.
### Calculate the mean intensity for each sample by chromosome.
inten.by.chrom <- meanIntensityByScanChrom(qxyData)
### Identify sex mis-annotation or sex chromosome aneuploidies.
mninten <- inten.by.chrom[[1]]
xcol <- rep(NA, nrow(scanAnnot))
xcol[scanAnnot$sex == "M"] <- "blue"
xcol[scanAnnot$sex == "F"] <- "red"
nx <- sum(snpAnnot$chromosome == autosomes+1) 
ny <- sum(snpAnnot$chromosome == autosomes+3) 
### Intensities
x1 <-mninten[,"X"]; y1 <- mninten[,"Y"]
main1 <- "Mean X vs \nMean Y Chromosome Intensity"
### Heterozygosity on X vs X intensity
x2 <- mninten[,"X"]; y2 <- scanAnnot$het.X
main2 <- "Mean X Chromosome Intensity vs
Mean X Chromosome Heterozygosity"
### Heterozygosity on X vs Y intensity
y3 <- mninten[,"Y"]; x3 <- scanAnnot$het.X
main3 <- "Mean X Chromosome Heterozygosity vs
Mean Y Chromosome Intensity"
### X vs A Heterozygosity
x4 <- scanAnnot$het.A[scanAnnot$sex == "F"]
y4 <- scanAnnot$het.X[scanAnnot$sex == "F"]
main4 <- "Mean Autosomal Heterozygosity vs
Mean X Chromosome Heterozygosity"
cols <- c("blue","red")
mf <- c("male", "female")
xintenlab <- paste("X intensity (n=", nx, ")", sep="")
yintenlab <- paste("Y intensity (n=", ny, ")", sep="")
pdf(paste (Misannotated_Sex_Check, "/01. Sex annotation analysis.pdf", sep=""), 7, 5)
par(mfrow=c(2,2))
plot(x1, y1, xlab=xintenlab, ylab=yintenlab,
     main=main1, col=xcol, cex.main=0.8)
legend("topright",mf,col=cols,pch=c(1,1))
plot(x2, y2, col=xcol, xlab=xintenlab, 
     ylab="X heterozygosity", main=main2, cex.main=0.8)
plot(x3, y3, col=xcol, ylab=yintenlab,
     xlab="X heterozygosity", main=main3, cex.main=0.8)
plot(x4,y4, col="red", xlab="Autosomal heterozygosity",
     ylab="X heterozygosity", main=main4, cex.main=0.8)
dev.off()


# Batch Quality Checks 
## Folder of for Batch Quality Checks outputs
path_Batch_Quality_Checks <- paste(path_quality_control,"/Batch_Quality_Checks", sep="") 
dir.create(path_Batch_Quality_Checks)
Missing_call_rate_before_filter <- paste(path_Batch_Quality_Checks,"/Missing_call_rate_before_filter", sep="") 
dir.create(Missing_call_rate_before_filter)
Missing_call_rate_after_filter <- paste(path_Batch_Quality_Checks,"/Missing_call_rate_after_filter", sep="") 
dir.create(Missing_call_rate_after_filter)

## Calculate number of missing calls per sex per snp over all samples
miss <- missingGenotypeBySnpSex(genoData)
#### Add missing.n1 to snpAnnot object
snpAnnot$missing.n1 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples",
  "except that females are excluded for Y chr SNPs")
#### Results:
##### 1. Missing counts by SNP				
###### Column with snpID
missing.counts <- data.frame(miss$missing.counts)
missing.counts <- data.frame (snpID=as.integer(rownames(missing.counts)), missing.counts)
###### Merge with snpName
missing.counts <- merge (x = snp_annot_data_frame_snpID[1:2], y =missing.counts, by = "snpID") 
###### Report
write.table(missing.counts, file = paste (Missing_call_rate_before_filter, 
                                          "/01. Missing counts by SNP and sex.txt", sep=""), sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 2. Sex count
write.table(data.frame(miss$scans.per.sex), file = paste (Missing_call_rate_before_filter,
                                                          "/02. Samples by sex.txt", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, 
            col.names = FALSE, quote = FALSE)
##### 3. Fraction of missing calls per SNP
###### Column with snpID
missing.fraction <- data.frame(missing.fraction = miss$missing.fraction)
missing.fraction <- data.frame (snpID=as.integer(rownames(missing.fraction)), missing.fraction)
###### Merge with snpName
missing.fraction <- merge (x = snp_annot_data_frame_snpID[1:2], y =missing.fraction, by = "snpID") 
###### Export the result
write.table(format(missing.fraction, digits = 2), file = paste (Missing_call_rate_before_filter,
                                                                "/03. Fraction of missing calls per SNP.txt", sep=""),
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
##### 4. Proportion of SNPs with MCR smaller than cutoff.
###### Function.
func_prop.SNPs.small.MCR <- function(chromosome, cutoof) {
  length(snpAnnot$missing.n1[snpAnnot$chromosome==chromosome]
         [snpAnnot$missing.n1[snpAnnot$chromosome==chromosome] <= cutoof]) /
    length(snpAnnot$missing.n1[snpAnnot$chromosome==chromosome])
}
###### Applying the function.
proportion_of_SNPs_above_missing_call_rate <- 
  sapply(seq( from=1, to = length(unique(snpAnnot$chromosome))+1), 
         func_prop.SNPs.small.MCR, cutoof=cutoof)
###### Merge proportion with chromosome
proportion_of_SNPs_above_missing_call_rate <-data.frame(
  "Chromosome" = seq( from=1, to = length(unique(snpAnnot$chromosome))+1),
  "Prop of SNPs above MCR" = proportion_of_SNPs_above_missing_call_rate)
###### Export the result
write.table(format(proportion_of_SNPs_above_missing_call_rate, digits = 2), 
            file = paste (Missing_call_rate_before_filter,
                          "/04. Proportion of SNPs above MCR by chromosome.txt", sep=""), sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 5. Number of SNPs in funciton of missing call rate
set.seed(1)
pdf(paste (Missing_call_rate_before_filter,
           "/05. Number of SNPs in function of MCR.pdf", sep=""), 7, 5)
hist(snpAnnot$missing.n1, ylim=c(0,num_SNP_auto), xlab="Missing call rate", 
     ylab = "Number of SNPs", main=" Number of SNPs in function of MCR")
dev.off()

### Calculate missing call rate per sample across SNPs.
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1] #exclude SNPs with 100% MCR
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
#### Add this vector to the annotation file
scanAnnot$missing.e1 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n1<1",
  "except that Y chr SNPs are excluded for females")
#### Calculate and store in the annotation file, the MCR for autosomes and the X chromosome.
auto <- colnames(miss$missing.counts) %in% 1:autosomes
missa <- rowSums(miss$missing.counts[,auto]) / sum(miss$snps.per.chr[auto])
missx <- miss$missing.counts[,"X"] / miss$snps.per.chr["X"]
scanAnnot$miss.e1.auto <- missa
scanAnnot$miss.e1.xchr <- missx
####  Create and store an object with identification of duplicated scans.
scanAnnot <- scanAnnot[order(scanAnnot$subjectID, scanAnnot$missing.e1),]
scanAnnot$duplicated <- duplicated(scanAnnot$subjectID)
##### Order scanAnnot by scanID
scanAnnot <- scanAnnot[order(scanAnnot$scanID),]
varMetadata(scanAnnot)["duplicated", "labelDescription"] <- "TRUE for duplicate scan with higher missing.e1"
#### Results: 
##### 6. Missing counts per sample by chromosome.
###### Create a column with scanID
missing.counts.sample.chromosome <- data.frame (miss$missing.counts)
missing.counts.sample.chromosome <- data.frame (
  scanID=as.integer(rownames(missing.counts.sample.chromosome)),
  missing.counts.sample.chromosome)
###### Merge with subjectID
missing.counts.sample.chromosome <- merge (
  x = Scan_Annotation_Data_Frame_scanID[1:2],
  y =missing.counts.sample.chromosome, by = "scanID") 
###### Assign column names
colnames(missing.counts.sample.chromosome) <- c("scanID", "subjectID", "Unknown", 
                                                paste("Chrom",seq(from=1, to =autosomes), sep = ""),
                                                "XChrom", "YChrom",
                                                "Mitochondrial")
###### Write the report
write.table(missing.counts.sample.chromosome, 
            file = paste (Missing_call_rate_before_filter,
                          "/06. Missing counts per sample by chromosome.txt", sep=""), sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 7. Missing SNPs per chromosome in all samples
missing.counts.chromosome <- data.frame (miss$snps.per.chr)
missing.counts.chromosome <- data.frame (chromosome=rownames(missing.counts.chromosome), 
                                         missing.counts.chromosome)
###### Writing the report
write.table(missing.counts.chromosome, 
            file = paste (Missing_call_rate_before_filter,
                          "/07. Missing SNPs per chromosome.txt", sep=""), sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 8. Missing fraction per sample
missing.fraction <- data.frame (miss$missing.fraction)
missing.fraction <- data.frame (
  scanID=as.integer(rownames(missing.fraction)), missing.fraction)
###### Merge with subjectID
missing.fraction <- merge (x = Scan_Annotation_Data_Frame_scanID[1:2],
                           y =missing.fraction, by = "scanID") 
###### Writing the report
write.table(format(missing.fraction, digits=2), 
            file = paste (Missing_call_rate_before_filter,
                          "/08. Missing fraction per sample.txt", sep=""), sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 9. Number of samples in funciton of missing call rate
set.seed(1)
pdf(paste (Missing_call_rate_before_filter,
           "/09. Number of samples in funciton of missing call rate.pdf", sep=""), 7, 5)
hist(scanAnnot$missing.e1, xlab="Missing call rate",
     ylab = "Number of samples",
     main="Number of samples in function of MCR")
dev.off()

### Calculate the MCR of the SNPs over samples whose MCR is greater than the cutoff of 0.05.
#### Vector of samples identifiers to remove
scan.exclude <- scanAnnot$scanID[scanAnnot$missing.e1 > 0.05]
#### Remove the samples
miss <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)
#### Calculate MCR
snpAnnot$missing.n2 <- miss$missing.fraction
#### Add metadata
varMetadata(snpAnnot)["missing.n2", "labelDescription"] <- 
  paste("fraction of genotype calls missing over all samples with missing.e1<0.05",
        "except that females are excluded for Y chr SNPs")

#### Results
##### 1. Missing counts by SNP				
###### Column with snpID
missing.counts <- data.frame(miss$missing.counts)
missing.counts <- data.frame (snpID=as.integer(rownames(missing.counts)), missing.counts)
###### Merge with snpName
missing.counts <- merge (x = snp_annot_data_frame_snpID[1:2], y =missing.counts, by = "snpID") 
###### Report
write.table(missing.counts, file = paste (Missing_call_rate_after_filter, 
                                          "/01. Missing counts by SNP and sex.txt", sep=""), sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)

##### 2. Sex count
write.table(data.frame(miss$scans.per.sex), file = paste (Missing_call_rate_after_filter,
                                                          "/02. Samples by sex.txt", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, 
            col.names = FALSE, quote = FALSE)
##### 3. Fraction of missing calls per SNP
###### Column with snpID
missing.fraction <- data.frame(missing.fraction = miss$missing.fraction)
missing.fraction <- data.frame (snpID=as.integer(rownames(missing.fraction)), missing.fraction)
###### Merge with snpName
missing.fraction <- merge (x = snp_annot_data_frame_snpID[1:2], y =missing.fraction, by = "snpID") 
###### Delimit number of decimal places
missing.fraction$missing.fraction <- format(missing.fraction$missing.fraction, digits = 2)
###### Export the result
write.table(format(missing.fraction, digits = 2), file = paste (Missing_call_rate_after_filter,
                                                                "/03. Fraction of missing calls per SNP.txt", sep=""),
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
##### 4. Proportion of SNPs with MCR smaller than cutoff.
###### Function.
func_prop.SNPs.small.MCR <- function(chromosome, cutoof) {
  length(snpAnnot$missing.n2[snpAnnot$chromosome==chromosome]
         [snpAnnot$missing.n2[snpAnnot$chromosome==chromosome] <= cutoof]) /
    length(snpAnnot$missing.n2[snpAnnot$chromosome==chromosome])
}
###### Applying the function.
proportion_of_SNPs_above_missing_call_rate <- 
  sapply(seq( from=1, to = length(unique(snpAnnot$chromosome))+1), 
         func_prop.SNPs.small.MCR, cutoof=cutoof)
###### Merge proportion with chromosome
proportion_of_SNPs_above_missing_call_rate <-data.frame(
  "Chromosome" = seq( from=1, to = length(unique(snpAnnot$chromosome))+1),
  "Prop of SNPs above MCR" = proportion_of_SNPs_above_missing_call_rate)
###### Export the result
write.table(format(proportion_of_SNPs_above_missing_call_rate, digits = 2), 
            file = paste (Missing_call_rate_after_filter,
                          "/04. Proportion of SNPs above MCR by chromosome.txt", sep=""), sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 5. Print and store a PDF whit the histogram of Missing Call Rate for All Probes
set.seed(1)
pdf(paste (Missing_call_rate_after_filter,
           "/05. Number of SNPs in function of MCR.pdf", sep=""), 7, 5)
hist(snpAnnot$missing.n2, ylim=c(0,num_SNP_auto), xlab="Missing call rate", 
     ylab = "Number of SNPs", main=" Number of SNPs in function of MCR")
dev.off()

### missing.e2
#### Vector of the SNPs to exclude.
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n2 >= 0.05]
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
#### Add the missing call rates vector to the sample annotation file
scanAnnot$missing.e2 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n2<0.05",
  "except that Y chr SNPs are excluded for females")
#### Results: 
##### 6. Missing counts per sample by chromosome.
###### Create a column with scanID
missing.counts.sample.chromosome <- data.frame (miss$missing.counts)
missing.counts.sample.chromosome <- data.frame (
  scanID=as.integer(rownames(missing.counts.sample.chromosome)),
  missing.counts.sample.chromosome)
###### Merge with subjectID
missing.counts.sample.chromosome <- merge (
  x = Scan_Annotation_Data_Frame_scanID[1:2],
  y =missing.counts.sample.chromosome, by = "scanID") 
###### Assign column names
colnames(missing.counts.sample.chromosome) <- c("scanID", "subjectID", "Unknown", 
                                                paste("Chrom",seq(from=1, to =autosomes), sep = ""),
                                                "XChrom", "YChrom",
                                                "Mitochondrial")
###### Write the report
write.table(missing.counts.sample.chromosome, 
            file = paste (Missing_call_rate_after_filter,
                          "/06. Missing counts per sample by chromosome.txt", sep=""), sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 7. Missing SNPs per chromosome in all samples
missing.counts.chromosome <- data.frame (miss$snps.per.chr)
missing.counts.chromosome <- data.frame (chromosome=rownames(missing.counts.chromosome), 
                                         missing.counts.chromosome)
###### Writing the report
write.table(missing.counts.chromosome, 
            file = paste (Missing_call_rate_after_filter,
                          "/07. Missing SNPs per chromosome.txt", sep=""), sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 8. Missing fraction per sample
missing.fraction <- data.frame (miss$missing.fraction)
missing.fraction <- data.frame (
  scanID=as.integer(rownames(missing.fraction)), missing.fraction)
###### Merge with subjectID
missing.fraction <- merge (x = Scan_Annotation_Data_Frame_scanID[1:2],
                           y =missing.fraction, by = "scanID") 
###### Writing the report
write.table(format(missing.fraction, digits=2), 
            file = paste (Missing_call_rate_after_filter,
                          "/08. Missing fraction per sample.txt", sep=""), sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
##### 9. Number of samples in funciton of missing call rate
set.seed(1)
pdf(paste (Missing_call_rate_after_filter,
           "/09. Number of samples in funciton of missing call rate.pdf", sep=""), 7, 5)
hist(scanAnnot$missing.e2, xlab="Missing call rate",
     ylab = "Number of samples",
     main="Number of samples in function of MCR")
dev.off()

## 1. Number of samples per batch.
set.seed(1)
pdf(paste (path_Batch_Quality_Checks, "/01. Number_of_samples_per_batch.pdf", sep=""), 7, 5)
barplot(table(scanAnnot$batch), ylab="Number of Samples", xlab="Batch", 
        main="Distribution of Samples per Batch")
dev.off()

## 2. MCR per batch in function of the number of samples of the batch
### Batches names
batches <- unique(scanAnnot$batch)
### Create an element to store the mean MCR per batch
bmiss <- rep(NA,length(batches)); names(bmiss) <- batches
### Create an element to store the length MCR per batch
bn <- rep(NA,length(batches)); names(bn) <- batches
### Calculate the means
for(i in 1:length(batches)) {
  mean_MCR_batch <- scanAnnot$missing.e1[is.element(scanAnnot$batch, batches[i])]
  bmiss[i] <- mean(mean_MCR_batch)
  bn[i] <- length(mean_MCR_batch)
}
### Generate a linear model to find the regression line of mean MCR per batch on the number of samples. 
lm_meanMCR_batch <- lm(bmiss ~ bn)
### Export plot
set.seed(1)
pdf(paste (path_Batch_Quality_Checks, "/02. Association between mean MCR per batch and number of samples.pdf",
           sep=""), 7, 5)
plot(bn, bmiss, xlab="Number of samples per batch", ylab="Mean MCR", 
     main="Mean Missing Call Rate vs\nSamples per Batch")
abline(lm_meanMCR_batch$coefficients)
dev.off()

## 3. Mean missing call rate per batch
write.table(data.frame(batch = unique(scanAnnot$batch), mean_MCR = format(mean_MCR_batch, digits = 2)), 
            file = paste (path_Batch_Quality_Checks, "/03. Mean MCR per batch", sep=""),
            sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)

## 4.	Chi-Square test of allelic frequency differences in batche
### 4.1. Association between batches and genetic origin
res <- batchChisqTest(genoData, batchVar="batch", return.by.snp=TRUE)
### Fraction of samples in each batch that are from each of the origins
batches <- unique(scanAnnot$batch)
eth_0 <- rep(NA,length(batches)); names(eth_0) <- sort(batches)
for(i in 1:length(batches)){
  x_0 <- scanAnnot$genetic_origin[is.element(scanAnnot$batch, batches[i])]
  xl_0 <- length(x_0[x_0 == 0]) # proportion of samples with genetic origin 0 in each batch
  eth_0[i] <- xl_0 / length(x_0) # mean of proportion
}
eth_1 <- rep(NA,length(batches)); names(eth_1) <- sort(batches)
for(i in 1:length(batches)){
  x_1 <- scanAnnot$genetic_origin[is.element(scanAnnot$batch, batches[i])]
  xl_1 <- length(x_1[x_1 == 1]) # proportion of samples with genetic origin 0 in each batch
  eth_1[i] <- xl_1 / length(x_1) # mean of proportion
}
set.seed(1)
pdf(paste (path_Batch_Quality_Checks, "/04.1. Association between batches and genetic origin.pdf",
           sep=""), 7, 5)
plot(eth_0,
     res$mean.chisq,
     xlab="Fraction of samples per batch of genetic origin 0",
     ylab="Average Chi-square Test Statistic",
     main="Fraction of samples per batch of 
     genetic origin 0 vs average chi-square test")
abline(v=mean(eth_0), lty=2, col="blue")
plot(eth_1,
     res$mean.chisq,
     xlab="Fraction of samples per batch of genetic origin 1",
     ylab="Average Chi-square Test Statistic",
     main="Fraction of samples per batch of 
     genetic origin 1 vs average chi-square test")
abline(v=mean(eth_1), lty=2, col="blue")
dev.off()
### 4.2. Mean chi-square and genomic inflation factor per batch 
write.table(data.frame(Batch = unique(scanAnnot$batch),
                       Mean.chi.square = format(res$mean.chisq, digits = 2),
                       Genomic.inflation.factor = format(res$lambda, digits = 2)),
            file = paste (path_Batch_Quality_Checks, 
                          "/04.2. Mean chi-square and genomic inflation factor per batch ", sep=""),
            sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)

#-------------------------------------------------------



#------------------ GWAS ---------------------

# The association analysis is performing excluding samples with low mean and/or median GC Score. 
# This is specified in the parameter scan.exclude.
GWAS<-assocRegression(genoData,
                      outcome = outcome,
                      model.type = model.type,
                      covar = covar,
                      ivar = ivar,
                      CI = CI,
                      scan.exclude = 
                        c(qual.results[qual.results$mean.quality < mean_GC.score_sample,][[1]]),
                      snpStart = 1,
                      snpEnd = num_SNP_auto,
                      block.size = block.size)
# Removing SNPs with low quality
## SNPs with low mean and/or median GC Score
GWAS_filt <- merge(x = data.frame(snpID = as.integer(as.character(qual.results_snp_filt$snpID))),
                   y = GWAS,
                   by = "snpID")
## SNPs associated with genetic origin
GWAS_filt <- subset(GWAS_filt, !(GWAS_filt$snpID %in% SNP_origin$snpID))
## SNPs with p-values in Hardy-Weinberg equilibrium test smaller than the cutoff
GWAS_filt <- merge(x = data.frame(snpID = hwe_filtered$snpID),
                   y = GWAS_filt,
                   by = "snpID")
## SNPs with a minor allele frequency smaller than the cutoff
GWAS_filt <- merge(x = data.frame(snpID = hwe_filtered_MAF$snpID),
                   y = GWAS_filt,
                   by = "snpID")
## SNPs with NA in p-val
GWAS_filt <- GWAS_filt[!is.na(GWAS_filt$Wald.pval)]

# Adjust P-values with multiple comparision
P_adjust <- p.adjust(GWAS_filt$Wald.pval, 
                     method = method,
                     n = length(GWAS$Wald.pval))
GWAS_P_adjust <- data.frame(GWAS_filt, P_adjust = P_adjust)

# Export results: 
## Create a folders for association analysis
association_analysis <- paste(path_outputs,"/association_analysis", sep="") 
dir.create(association_analysis)
## Export association dataframe
write.table(GWAS_P_adjust,
            file = paste (association_analysis, "/", "Association result.txt", sep=""), 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)
## Manhattan plot
pdf(paste (association_analysis, "/", "Manhattan_Plot.pdf", sep=""), 7, 5)
manhattanPlot(p = GWAS_P_adjust$P_adjust,
              chromosome = GWAS_P_adjust$chr,
              signif = signif)
dev.off()
## SNPs associated list
Associated_SNPs <- GWAS_P_adjust[GWAS_P_adjust$P_adjust <= signif,]
### Remove NAs
Associated_SNPs <- Associated_SNPs[!is.na(Associated_SNPs$Wald.pval),]
### Complete data frame with more information of the SNP
Associated_SNPs <- merge(x = Associated_SNPs, y = snp_annot_data_frame_snpID, by = "snpID")
write.table(Associated_SNPs, file = paste (association_analysis, "/", "Associated_SNPs.txt", sep=""),
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE,
            quote = FALSE)
