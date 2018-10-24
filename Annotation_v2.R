# Gene anotation

# Packages
library(mygene)

# Parameter
  # Gene ID in Entrez/ensembl
gene = 282659

# Anotation
Annotation <- getGene(gene, fields="all", return.as = "records")

# Export information

## Create the results's folder
Results <- "../Results" 
dir.create(Results)

## 1. Gene description
Gene_description <- data.frame (x= c ("Entrez/ensembl gene id:", "name:", "symbol:", "other_names:",
                                      "Type_of_gene:", "Taxonomy_id:", "ensembl type of gene:",
                                      "chromosome:", "start position:", "end_position:"),
                               y = c(Annotation[[1]][["_id"]], Annotation[[1]][["name"]],
                                     Annotation[[1]][["symbol"]], 
                                     paste(Annotation[[1]][["other_names"]], collapse = "; "),
                                     Annotation[[1]][["type_of_gene"]], Annotation[[1]][["taxid"]],
                                     Annotation[[1]][["ensembl"]][["type_of_gene"]],
                                     Annotation[[1]][["genomic_pos"]][["chr"]],
                                     Annotation[[1]][["genomic_pos"]][["start"]],
                                     Annotation[[1]][["genomic_pos"]][["end"]])
                               )
### Export Gene description
write.table(Gene_description, file = paste(Results, "/", "1. Gene_description.txt", sep = ""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

## 2. Data bases ids
Data_bases_Ids <- data.frame (
  data_base = c ("Entrez/ensembl gene id", "UniGene", "GenBank_genomic", "GenBank_protein",
                 "GenBank_rna", "ensembl_gene", "ensembl_protein", "ensembl_transcript", "Pfam",
                 "Prosite", "Uniprot"),
  Id = c(Annotation[[1]][["_id"]], Annotation[[1]][["unigene"]],
         paste(Annotation[[1]][["accession"]][["genomic"]], collapse = "; "),
         paste(Annotation[[1]][["accession"]][["protein"]], collapse = "; "),
         paste(Annotation[[1]][["accession"]][["rna"]], collapse = "; "),
         paste(Annotation[[1]][["ensembl"]][["gene"]], collapse = "; "),
         paste(Annotation[[1]][["ensembl"]][["protein"]], collapse = "; "),
         paste(Annotation[[1]][["ensembl"]][["transcript"]], collapse = "; "),
         paste(Annotation[[1]][["pfam"]], collapse = "; "),
         paste(Annotation[[1]][["prosite"]], collapse = "; "),
         paste(Annotation_multiple_genes[[2]][["uniprot"]], collapse = "; ")))

### Export Data_bases_Ids
write.csv(Data_bases_Ids, file = paste(Results, "/", "2. Data_bases_Ids.csv", sep = ""),
          row.names = FALSE)

## 3. Translations
### 3.1 translations GeneBank ids
#### vector with rna ids
rna=0
for (i in seq(from = 1, to = length(Annotation[[1]][["accession"]][["translation"]]))) {
  rna[i] = Annotation[[1]][["accession"]][["translation"]][[i]][[2]]
}
#### vector with protein ids
protein=0
for (i in seq(from = 1, to = length(Annotation[[1]][["accession"]][["translation"]]))) {
  protein[i] = Annotation[[1]][["accession"]][["translation"]][[i]][[1]]
}
#### Data frame
translations_GeneBank_ids <- data.frame(rna = rna,
                           protein = protein)
#### Export translations_GeneBank_ids
write.csv(translations_GeneBank_ids, file = paste(Results, "/", "3.1. translations_GeneBank_ids.csv", sep = ""),
          row.names = FALSE)

### 3.2 Translations ensembl ids
#### vector with rna ids
rna=0
for (i in seq(from = 1, to = length(Annotation[[1]][["ensembl"]][["translation"]]))) {
  rna[i] = Annotation[[1]][["accession"]][["translation"]][[i]][[2]]
}
#### vector with protein ids
protein=0
for (i in seq(from = 1, to = length(Annotation[[1]][["ensembl"]][["translation"]]))) {
  protein[i] = Annotation[[1]][["accession"]][["translation"]][[i]][[1]]
}
#### Data frame
translations_ensembl_ids <- data.frame(rna = rna,
                                        protein = protein)
#### 3.2 Export translations_ensembl_ids
write.csv(translations_ensembl_ids, file = paste(Results, "/", "3.2. translations_ensembl_ids.csv", sep = ""),
          row.names = FALSE)

## 4. Exons
### End of the coding region
cdsend = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  cdsend[i] = Annotation[[1]][["exons"]][[i]][["cdsend"]]
}
### Start of the coding region
cdsstart = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  cdsstart[i] = Annotation[[1]][["exons"]][[i]][["cdsstart"]]
}
### Chromosome
chr = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  chr[i] = Annotation[[1]][["exons"]][[i]][["chr"]]
}
### Strand
strand = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  strand[i] = Annotation[[1]][["exons"]][[i]][["strand"]]
}
### Transcript
transcript = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  transcript[i] = Annotation[[1]][["exons"]][[i]][["transcript"]]
}
### Transcription end position
txend = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  txend[i] = Annotation[[1]][["exons"]][[i]][["txend"]]
}
### Transcription start position
txstart = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["exons"]]))) {
  txstart[i] = Annotation[[1]][["exons"]][[i]][["txstart"]]
}
### Data frame
exons <- data.frame(Transcript_Id = transcript, Chromosome = chr, Strand = strand, 
                    Transcription_start_position = txstart, Transcription_end_position = txend,
                    Coding_region_start = cdsstart, Coding_region_end = cdsend)
### Export exons
write.csv(exons, file = paste(Results, "/", "4. exons.csv", sep = ""),
          row.names = FALSE)

## 5. Pubmed
### Pubmed Id
pmid = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["generif"]]))) {
  pmid[i] = Annotation[[1]][["generif"]][[i]][["pubmed"]]
}
### Description
text = 0
for (i in seq(from = 1, to = length(Annotation[[1]][["generif"]]))) {
  text[i] = Annotation[[1]][["generif"]][[i]][["text"]]
}
### Data frame
pubmed <- data.frame(Pubmed_Id = pmid,
                     Description = text)
### Export pubmed
write.csv(pubmed, file = paste(Results, "/", "5. pubmed.csv", sep = ""),
          row.names = FALSE)

## 6. Gene ontology
### Biological Process (BP)
#### Initialize variables: Experimental evidence, id, and term
evidence_bp = 0
id_bp = 0
term_bp = 0
#### Asign vectors to variables. The function identify if there are more than one BP
if (length(Annotation[[1]][["go"]][["BP"]][[1]]) == 1) {
  evidence_bp = Annotation[[1]][["go"]][["BP"]][["evidence"]]
  id_bp = Annotation[[1]][["go"]][["BP"]][["id"]] 
  term_bp = Annotation[[1]][["go"]][["BP"]][["term"]] } else
    for (i in seq(from = 1, to = length(Annotation[[1]][["go"]][["BP"]]))) {
      evidence_bp[i] = Annotation[[1]][["go"]][["BP"]][[i]][["evidence"]] 
      id_bp[i] = Annotation[[1]][["go"]][["BP"]][[i]][["id"]]
      term_bp[i] = Annotation[[1]][["go"]][["BP"]][[i]][["term"]]
    }
### Cellular Component (CC)
#### Initialize variables: Experimental evidence, id, and term
evidence_cc = 0
id_cc = 0
term_cc = 0
#### Asign vectors to variables. The function identify if there are more than one CC
if (length(Annotation[[1]][["go"]][["CC"]][[1]]) == 1) {
  evidence_cc = Annotation[[1]][["go"]][["CC"]][["evidence"]]
  id_cc = Annotation[[1]][["go"]][["CC"]][["id"]] 
  term_cc = Annotation[[1]][["go"]][["CC"]][["term"]] } else
    for (i in seq(from = 1, to = length(Annotation[[1]][["go"]][["CC"]]))) {
      evidence_cc[i] = Annotation[[1]][["go"]][["CC"]][[i]][["evidence"]] 
      id_cc[i] = Annotation[[1]][["go"]][["CC"]][[i]][["id"]]
      term_cc[i] = Annotation[[1]][["go"]][["CC"]][[i]][["term"]]
    }
### Molecular Function (MF)
#### Initialize variables: Experimental evidence, id, and term
evidence_mf = 0
id_mf = 0
term_mf = 0
#### Asign vectors to variables. The function identify if there are more than one MF
if (length(Annotation[[1]][["go"]][["MF"]][[1]]) == 1) {
  evidence_mf = Annotation[[1]][["go"]][["MF"]][["evidence"]]
  id_mf = Annotation[[1]][["go"]][["MF"]][["id"]] 
  term_mf = Annotation[[1]][["go"]][["MF"]][["term"]] } else
    for (i in seq(from = 1, to = length(Annotation[[1]][["go"]][["MF"]]))) {
      evidence_mf[i] = Annotation[[1]][["go"]][["MF"]][[i]][["evidence"]] 
      id_mf[i] = Annotation[[1]][["go"]][["MF"]][[i]][["id"]]
      term_mf[i] = Annotation[[1]][["go"]][["MF"]][[i]][["term"]]
    }
### Data frame
GO <- data.frame(category = c(rep("Biological Process",length(Annotation[[1]][["go"]][["BP"]][[1]])),
                              rep("Cellular Component",length(Annotation[[1]][["go"]][["CC"]][[1]])),
                              rep("Molecular Function",length(Annotation[[1]][["go"]][["MF"]][[1]]))),
                 experimental_evidence = c(evidence_bp, evidence_cc, evidence_mf),
                 id = c(id_bp, id_cc, id_mf),
                 term = c(term_bp, term_cc, term_mf))
### Export 6. GO
write.csv(GO, file = paste(Results, "/", "6. Gene Ontology.csv", sep = ""),
          row.names = FALSE)
