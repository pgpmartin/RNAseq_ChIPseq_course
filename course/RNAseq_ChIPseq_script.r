## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)

## basedir="/work/pmartin/PROJECTS/RNAseqChIPseq"

## cd ${basedir}/bank

## cd ${basedir}/bank

## cd ${basedir}/bank

## awk '$1 ~ /^##.+/ || $1=="chr19"|| $1=="chr20"' \

## mkdir -p ${basedir}/bank/star_index

## mkdir -p ${basedir}/bank/bowtie2_index

## cd ${basedir}/bank

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/008/SRR2084598/SRR2084598.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR208/008/SRR2084598

## ---- eval=FALSE---------------------------------------------------------
## #All the files:
## gse70408 <- read.table(file.path("../data","PRJNA288695.txt"),
##                        sep="\t", header=T)
## 
## #Selected files:
## selgse70408 <- gse70408[c(
##                     grep("total.RNA", gse70408[,2]),
##                     grep("PAF1_WT", gse70408[,2])),]
## 
## #Format the links for aspera client:
## selgse70408$fasp_link <- gsub("ftp.sra.ebi.ac.uk",
##                               "era-fasp@fasp.sra.ebi.ac.uk:",
##                               selgse70408$fastq_ftp)
## 
## #save the links
## write.table(selgse70408$fasp_link,
##             file.path("../data/downloadSelectedGSE70408.txt"),
##             row.names = FALSE, col.names = FALSE, quote = FALSE)

## cat downloadSelectedGSE70408.txt | while read LIST

## mkdir -p ${basedir}/RNAseq/rawQC

## mkdir -p ${basedir}/RNAseq/Aligned

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
## #module load system/R-3.4.3 #on genologin cluster (choose the right version for you)
## library(GenomicAlignments)
## #Prepare the annotations:
## library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## txdb <- GenomeInfoDb::keepSeqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene,
##                                     c("chr19", "chr20"),
##                                     pruning.mode = "coarse")
## hg38exBygn <- exonsBy(txdb, by="gene")
## 
## # Bam files:
## bampath <- "/work/pmartin/PROJECTS/RNAseqChIPseq/RNAseq/Aligned"
## bamfn <- BamFileList(dir(bampath, pattern="*.sortedByCoord.out.bam$", full=TRUE))
## 
## ct <- summarizeOverlaps(hg38exBygn,
##                         bamfn,
##                         mode = "Union",
##                         singleEnd=TRUE,
##                         ignore.strand = TRUE)
## colnames(ct) <- paste0("SRR208459",6:9)

## ---- eval=FALSE, include=FALSE------------------------------------------
## saveRDS(head(assays(ct)$counts), "../RData/headct.rds")
## saveRDS(ct, "../RData/ct.rds")

## ---- eval=FALSE---------------------------------------------------------
## head(assays(ct)$counts)

## ---- echo=FALSE---------------------------------------------------------
headct <- readRDS("../RData/headct.rds")
headct %>%
    kable() %>%
    kable_styling(font_size = 13,
                  bootstrap_options = c("striped", "condensed"),
                  full_width = FALSE)

## ---- eval=FALSE---------------------------------------------------------
## 100 * colSums(assays(ct)$counts) / countBam(bamfn)$records
## # ~76-80%

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
## library(DESeq2)
## colData <- DataFrame(SampleName = colnames(ct),
##                      condition = factor(rep(c("shPAF1","Ctrol"), each=2)),
##                      row.names = colnames(ct))
## 
## dds <- DESeqDataSetFromMatrix(countData = assays(ct)$counts,
##                               colData = colData,
##                               design = ~ condition)

## ---- eval=FALSE---------------------------------------------------------
## dds <- dds[rowSums(assay(dds))!=0,]

## ---- eval=FALSE, message=FALSE------------------------------------------
## dds <- DESeq(dds)

## ---- eval=FALSE, include=FALSE------------------------------------------
## saveRDS(results(dds), "../RData/ddsres.rds")

## ---- eval=FALSE---------------------------------------------------------
## sum(results(dds)$padj<0.05, na.rm = TRUE)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
ddsres <- readRDS("../RData/ddsres.rds")
sum(ddsres$padj<0.05, na.rm=TRUE)

## ---- eval=FALSE---------------------------------------------------------
## head(results(dds))

## ---- echo=FALSE---------------------------------------------------------
head(ddsres) %>%
    kable(digits=c(2,2,3,3,6,6)) %>%
    kable_styling(font_size = 13,
                  bootstrap_options = c("striped", "condensed"),
                  full_width = FALSE)

## ---- eval=FALSE, include=FALSE------------------------------------------
## ## Compare summarizeOverlaps with STAR counting:
## gnct <- lapply(dir(pattern="_ReadsPerGene.out.tab"), read.table, sep="\t")
## ctstar <- do.call("cbind", lapply(gnct, `[`, -(1:4), 2))
##   rownames(ctstar) <- gsub("\\..+", "", gnct[[1]][-(1:4),1])
##   colnames(ctstar) <- paste0("SRR208459",6:9)
## 
##   #Problem, the names don't match
## library(org.Hs.eg.db)
## eg2ens <- select(org.Hs.eg.db, keytype="ENTREZID", column="ENSEMBL", key=rownames(ct))
## #There are 40 Entrez gene that match with 2 ensembl genes (also because Ensembl database version is different between star counting and org package)
## #let's keep only the single match (also for Ensembl):
## egdup <- eg2ens[duplicated(eg2ens[,"ENTREZID"]),"ENTREZID"]
## eg2ens <- eg2ens[!(eg2ens[,"ENTREZID"] %in% egdup),]
## ensdup <- eg2ens[duplicated(eg2ens[,"ENSEMBL"]),"ENSEMBL"]
## eg2ens <- eg2ens[!(eg2ens[,"ENSEMBL"] %in% ensdup), ]
## eg2ens <- eg2ens[eg2ens[,"ENSEMBL"] %in% rownames(ctstar), ]
## 
## #get non ambiguous count tables
## ctstarOK <- ctstar[eg2ens[,"ENSEMBL"],]
## ctOK <- ct[eg2ens[,"ENTREZID"],]
## 
## par(mfrow=c(2,2))
## for (i in 1:4) { plot(log2(ctstarOK[,i]+1), log2(assays(ctOK)$count[,i]+1)) }
## #OK counts are really close. I should use the exact same annotation to do this comparison
## 
## ##We should use these annotations:
## basedir="/work/pmartin/PROJECTS/RNAseqChIPseq"
## txdb <- makeTxDbFromGFF(file.path(basedir, "bank", "gencode.v28.annotation_chr1920.gtf"))

## ----Estimate gene length, eval=FALSE------------------------------------
## txsByGene <- transcriptsBy(txdb, "gene")
## GeneLength <- median(width(txsByGene))
## GeneLengthInKb <- GeneLength / 1e3

## ----RNAseq RPKM values, eval=FALSE--------------------------------------
## #geneLengths is in Kb
## rpkm_fun <- function(countMatrix, geneLengths) {
##     x <- t( t(countMatrix) * 1e6 / colSums(countMatrix) ) # get RPM
##     rpkm.mat <- x / geneLengths # scale for gene length
##     return(rpkm.mat)
## }
## 
## RPKM <- rpkm_fun(assays(ct)$counts,
##                  GeneLengthInKb)

## ----RNAseq RPKM Ryan Devon function, eval=FALSE, include=FALSE----------
## #The function above is strictly identical to this function proposed by Ryan Devon (https://www.biostars.org/p/96084/):
## rpkm_fun <- function(countMatrix, geneLengths) {
##     (10^6)*t(t(countMatrix/geneLengths)/colSums(countMatrix))
## }
## #But better illustrates the 2 steps that are inversed in TPM

## ----saveRPKM, eval=FALSE, include=FALSE---------------------------------
## saveRDS(RPKM, "../RData/RPKM.rds")

## ----edgeR RPKM, eval=FALSE----------------------------------------------
## edg <- DGEList(assays(ct)$counts,
##                group = factor(rep(c("shPAF1","Ctrol"), each=2)))
## edg <- edgeR::calcNormFactors(edg)
## 
## edgRPKM <- rpkm(edg,
##                 gene.length = GeneLength,
##                 normalized.lib.sizes = FALSE,
##                 log = FALSE,
##                 prior.count = 0)
## 
## identical(round(RPKM,10), round(edgRPKM, 10)) #TRUE

## ----RNAseq TPM values, eval=FALSE---------------------------------------
## #geneLengths is in Kb
## tpm_fun <- function(countMatrix, geneLengths) {
##     x <- countMatrix / geneLengths #get RPK
##     tpm.mat <- t( t(x) * 1e6 / colSums(x) ) #Scale for sequencing depth
##     return(tpm.mat)
## }
## 
## TPM <- tpm_fun(assays(ct)$counts,
##                GeneLengthInKb)

## ----saveTPM, eval=FALSE, include=FALSE----------------------------------
## saveRDS(TPM, "../RData/TPM.rds")

## ---- eval=FALSE---------------------------------------------------------
## mybam <- readGAlignments(file.path(bampath,
##                                    "SRR2084596_chr1920_Aligned.sortedByCoord.out.bam"))
## # !! for RNA-seq do not convert the reads to granges because spliced reads would results in coverage on the introns
## covbam <- coverage(grglist(mybam))

## ----RNAseq export profile, eval=FALSE, message=FALSE, warning=FALSE-----
## rtracklayer::export(covbam, "myProfile.bigwig")

## zcat ${basedir}/data/SRR208453{7..9}.fastq.gz > ${basedir}/data/PAF1_ChIPseq.fastq

## mkdir -p ${basedir}/ChIPseq/rawQC

## mkdir -p ${basedir}/ChIPseq/Aligned

## ----ChIPseq sam2Sortedbam, eval=FALSE-----------------------------------
## module load bioinfo/samtools-1.4 #on genologin
## samtools view -b ${outfile} | samtools sort --threads 4 -o ${outfile%%.sam}.bam -

## module load bioinfo/picard-2.14.1 #sets the ${PICARD} variable on genologin

## ----ChIPseq index aligned read, eval=FALSE------------------------------
## AlignedReads="${outfile%%.sam}_MarkedDup.bam"
## samtools index ${AlignedReads}

## ----ChIPseq stats on aligned reads, eval=FALSE--------------------------
## #aligned reads by chromosome using idxstats
## samtools idxstats ${AlignedReads} > ${AlignedReads%%.bam}.idxstats
## 
## #QC passed/failed reads using flagstat
## samtools flagstat ${AlignedReads} > ${AlignedReads%%.bam}.flagstat
## 
## #More statistics (here removing duplicates for illustrating purposes)
## samtools stats \
##   --remove-dups \
## ${AlignedReads} > \
## ${AlignedReads%%.bam}.bamstats
## 
## #Plots  from these statistics can be generated using the plot-bamstats function present in the misc subfolder of samtools

## ----ChIPseq coverage in R, eval=FALSE, message=FALSE, warning=FALSE-----
## # cd ${basedir}/ChIPseq/Aligned
## # module load system/R-3.4.3 #on genologin cluster (choose the right version for you)
## library(chipseq)
## library(GenomicAlignments)
## 
## #import BAM
## paf1bam <- readGAlignments("PAF1_ChIPseq_chr1920_MarkedDup.bam")
## #estimate fragment length
## fraglen <- chipseq::estimate.mean.fraglen(granges(paf1bam)) #method="correlation" for SPP-like
## #get coverage on extended reads (~200bp)
## paf1cov <- coverage(resize(granges(paf1bam), 200))

## ----export PAF1 coverage, eval=FALSE, message=FALSE, warning=FALSE------
## library(rtracklayer)
## rtracklayer::export(paf1cov, "PAF1_coverage.bigwig")

## ----saveRDS paf1cov, eval=FALSE-----------------------------------------
## saveRDS(paf1cov, "../RData/paf1cov.rds")

## module load system/Python-2.7.2 #on genologin cluster

## ----import macs2 peaks in R, eval=FALSE, message=FALSE, warning=FALSE----
## # cd ${basedir}/ChIPseq/Aligned
## # module load system/R-3.4.3 #on genologin cluster (choose the right version for you)
## paf1Peaks <- genomation::readNarrowPeak("macs2_PAF1_peaks.narrowPeak")

## ----import gtf annotations in R, eval=FALSE, message=FALSE--------------
## basedir="/work/pmartin/PROJECTS/RNAseqChIPseq"
## txdb <- GenomicFeatures::makeTxDbFromGFF(file.path(basedir,
##                                                    "bank",
##                                                    "gencode.v28.annotation_chr1920.gtf"))

## ----ChIPseq PAF1 peak annotation ChIPseeker, eval=FALSE, message=FALSE, warning=FALSE----
## paf1peakanno <- ChIPseeker::annotatePeak(paf1Peaks,
##                                          tssRegion=c(-1000, 1000),
##                                          TxDb=txdb)

## ----saveRDS paf1peakanno, eval=FALSE, include=FALSE---------------------
## saveRDS(paf1peakanno, "../RData/paf1peakanno.rds")

## ----readRDS paf1peakanno, eval=FALSE, include=FALSE---------------------
## paf1peakanno <- readRDS("../RData/paf1peakanno.rds")

## ----ChIPseq plotAnnoBar, fig.height=3, fig.align = "center", message=FALSE, eval=FALSE----
## #png("course/figs/plotAnnoBar_paf1peakanno.png", width=750, height=450, res=150)
## ChIPseeker::plotAnnoBar(paf1peakanno)
## #dev.off()

## ----ChIPseq upsetplot, fig.height=5, fig.width=10, fig.align = "center", eval=FALSE----
## #png("course/figs/upset_paf1peakanno.png", width=1500, height=750, res=150)
## ChIPseeker::upsetplot(paf1peakanno, vennpie=TRUE)
## #dev.off()

## ----readRDS paf1cov, eval=FALSE, include=FALSE--------------------------
## paf1cov <- readRDS("../RData/paf1cov.rds")

## ----ChIPseq get TSS region, eval=FALSE----------------------------------
## txBygn <- GenomicFeatures::transcriptsBy(txdb, "gene")
## tss <- unique(unlist(promoters(txBygn, upstream=0, downstream=1)))

## ----ChIPseq get PAF1 profile around TSS, eval=FALSE---------------------
## paf1prof <- paf1cov[tss + 1000]
## names(paf1prof) <- tss$tx_id
## paf1prof[strand(tss)=="-"] <- lapply(paf1prof[strand(tss)=="-"], rev)

## ----saveRDS paf1prof, eval=FALSE, include=FALSE-------------------------
## saveRDS(paf1prof, "../RData/paf1prof.rds")

## ----readRDS paf1prof, eval=TRUE, include=FALSE--------------------------
paf1prof <- readRDS("../RData/paf1prof.rds")

## ----ChIPseq convert paf1prof to matrix, eval=FALSE----------------------
## paf1profmat <- matrix(as.numeric(unlist(paf1prof, use.names=FALSE)),
##                                  nrow = length(paf1prof),
##                                  byrow = TRUE,
##                                  dimnames = list(names(paf1prof), NULL))

## ----ChIPseq select genes for heatmap, eval=FALSE------------------------
## avgpaf1 <- apply(paf1profmat[,1001:1501], 1, mean)
## paf1profmat_sel <- paf1profmat[order(avgpaf1, decreasing = TRUE),][1:500,]

## ----ChIPseq average signal in bins of 10bp, eval=FALSE------------------
## paf1profmat_sel <- t(apply(paf1profmat_sel, 1, tapply,
##                          c(rep(1:100, each = 10),
##                            101,
##                            rep(102:201, each = 10)),
##                          mean))

## ----ChIPseq prepare TSS mat, eval=FALSE---------------------------------
## paf1TSSmat <- paf1profmat_sel[,-101]
## attr(paf1TSSmat, "upstream_index") <- 1:100
## attr(paf1TSSmat, "target_index") <- numeric(0)
## attr(paf1TSSmat, "downstream_index") <- 101:200
## attr(paf1TSSmat, "extend") <- c(1000,1000)
## attr(paf1TSSmat, "signal_name") <- "PAF1"
## attr(paf1TSSmat, "target_name") <- "TSS"
## attr(paf1TSSmat, "target_is_single_point") <- TRUE
## class(paf1TSSmat) = c("normalizedMatrix", "matrix")

## ----ChIPseq saveRDS paf1TSSmat, eval=FALSE, include=FALSE---------------
## saveRDS(paf1TSSmat, "../RData/paf1TSSmat.rds")

## ----ChIPseq readRDS paf1TSSmat, eval=FALSE, include=FALSE---------------
## paf1TSSmat <- readRDS("../RData/paf1TSSmat.rds")

## ----ChIPseq heatmap, message=FALSE, warning=FALSE, fig.width = 3, fig.align = "center", eval=FALSE----
## library(EnrichedHeatmap)
## col_fun = circlize::colorRamp2(c(0, 30), c("white", "red"))
## #png("course/figs/hmap_paf1TSSmat.png", width=450, height=750, res=150)
##     EnrichedHeatmap(paf1TSSmat,
##                 row_order=1:nrow(paf1TSSmat),
##                 col=col_fun,
##                 name="PAF1",
##                 axis_name=c('-1Kb','TSS','+1Kb'))
## #dev.off()

## ----ChIPseq calculate average profile, eval=FALSE-----------------------
## avgprof <- apply(paf1profmat, 2, mean)

## ----ChIPseq saveRDS avgprof, eval=FALSE, include=FALSE------------------
## saveRDS(avgprof, "../RData/avgprof.rds")

## ----ChIPseq readRDS avgprof, eval=TRUE, include=FALSE-------------------
avgprof <- readRDS("../RData/avgprof.rds")

## ----ChIPseq plot avgprof, eval=TRUE, fig.align = "center"---------------
plot(-1000:1000,
     avgprof, 
     type="l",
     lwd=3,
     axes=FALSE,
     xlab="Position relative to TSS",
     ylab="Average coverage (reads)")
axis(side=1,
     at=c(-1000, 0, 1000),
     labels = c("-1Kb", "TSS", "+1Kb"))
axis(side = 2)
abline(v=0, lty=2, lwd=1.5, col="darkgrey")

## ----CHIPseq select tssWithPAF1, eval=FALSE------------------------------
## tssWithPAF1 <- tss[overlapsAny(tss+1000, paf1Peaks)]

## ---- load BSgenome.Hsapiens.UCSC.hg38, eval=FALSE-----------------------
## library(BSgenome.Hsapiens.UCSC.hg38)

## ----ChIPseq get core promoters of tssWithPAF1, eval=FALSE---------------
## selseq <- BSgenome::getSeq(Hsapiens,
##                            promoters(tssWithPAF1,
##                                      upstream=500,
##                                      downstream=1))

## ----ChIPseq export tssWithPAF1 promoters to fasta, eval=FALSE-----------
## writeXStringSet(selSeq, "../data/Prom_tssWithPAF1.fa")

