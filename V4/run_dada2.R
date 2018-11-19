#!/usr/bin/Rscript

# The following files need to be in the work_dir folder:
# - silva_nr_v123_train_set.fa.gz
# - silva_species_assignment_v123.fa.gz

library(dada2)
library(stringr)
library(reshape2)

data_dir = "/media/harddrive/sander/DADA2_Eline_070617/data"
work_dir = "/media/harddrive/sander/DADA2_Eline_070617/"

setwd(work_dir)

# Making list of libraries and their files

print("Making list of libraries and their files")
print("")

files = list.files(data_dir)
pattern = "([a-zA-Z0-9\\-]+_[a-zA-Z0-9]+)_[a-zA-Z0-9]+_(R1|R2)_[0-9]+.fastq.gz"
files = str_match(files, pattern)
files = data.frame(files)
names(files) = c("fileName", "libraryName", "read")
files = files[!is.na(files$fileName),]
files = files[files$libraryName!="Undetermined_S0",]
files = dcast(files, formula=libraryName~read, value.var="fileName")

filenameFs = file.path(data_dir, files$R1)
filenameRs = file.path(data_dir, files$R2)

save(files, file="files.Robject")

# Paired read filtering and trimming

print("Paired read filtering and trimming")
print("")

filt_dir = file.path(data_dir, "filtered")
if (!file_test("-d", filt_dir)) dir.create(filt_dir)

filenameFiltFs = file.path(filt_dir, str_c(files$libraryName, "_F_filtered.fastq.gz"))
filenameFiltRs = file.path(filt_dir, str_c(files$libraryName, "_R_filtered.fastq.gz"))

for (i in 1:nrow(files)) {
  fastqPairedFilter(fn=c(filenameFs[i], filenameRs[i]),
                    fout=c(filenameFiltFs[i], filenameFiltRs[i]),
                    maxN=0, truncQ=2, maxEE=2,
                    trimLeft=c(12, 12), truncLen=c(240, 220),
                    verbose=T, compress=T)
}

# Some of the original data files could have been empty so we need to re-construct the filtered file names
filenameFiltFs = list.files(filt_dir, pattern="_F_filtered.fastq.gz", full.names=T)
filenameFiltRs = list.files(filt_dir, pattern="_R_filtered.fastq.gz", full.names=T)

# Dereplicating reads

print("Dereplicating reads")
print("")

derepFs = derepFastq(filenameFiltFs, verbose=F)
derepRs = derepFastq(filenameFiltRs, verbose=F)
names(derepFs) = str_match(filenameFiltFs, "\\/([a-zA-Z0-9_-]+)_F_filtered\\.fastq\\.gz")[,2]
names(derepRs) = str_match(filenameFiltRs, "\\/([a-zA-Z0-9_-]+)_R_filtered\\.fastq\\.gz")[,2]

save(derepFs, file="derepFs.Robject")
save(derepRs, file="derepRs.Robject")

# Applying dada error correction

print("Applying dada error correction")
print("")

dadaFs = dada(derepFs, selfConsist=T, err=NULL, multithread=T)
dadaRs = dada(derepRs, selfConsist=T, err=NULL, multithread=T)

save(dadaFs, file="dadaFs.Robject")
save(dadaRs, file="dadaRs.Robject")

# Merging reads and making sequence table

print("Merging reads and making sequence table")
print("")

mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)
save(mergers, file="mergers.Robject")
seqTable = makeSequenceTable(mergers)
print(table(str_length(colnames(seqTable))))

save(seqTable, file="seqTable.Robject")

# Removing too long sequences

print("Removing too long sequences")
print("")

seqTable = seqTable[,str_length(colnames(seqTable)) <= 230]
print(table(str_length(colnames(seqTable))))

# Removing chimera's

print("Removing chimera's")
print("")

seqTable.nochim = removeBimeraDenovo(seqTable, verbose=TRUE)

readCounts = colSums(seqTable.nochim)
readLengths = str_length(colnames(seqTable.nochim))
print(tapply(readCounts, readLengths, sum))

seqTable.nochim.fil = seqTable.nochim[,str_length(colnames(seqTable.nochim))<=230]

save(seqTable.nochim.fil, file="seqTable.nochim.fil.Robject")

# Assigning taxonomy

print("Assigning taxonomy")
print("")

taxonTable = assignTaxonomy(seqTable.nochim.fil, refFasta="silva_nr_v123_train_set.fa.gz")
taxonTable = addSpecies(taxonTable, refFasta="silva_species_assignment_v123.fa.gz", allowMultiple=T)
colnames(taxonTable) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

save(taxonTable, file="taxonTable.Robject")
