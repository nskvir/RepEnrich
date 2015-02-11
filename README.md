#RepEnrich
## Tutorial By Steven Criscione
Email: [steven_criscione@brown.edu](mailto:steven_criscione@brown.edu)

### Dependencies
This example is for mouse genome **mm9**. Before
getting started you should make sure you have installed the dependencies
for RepEnrich. RepEnrich requires python version 2.7.3.
RepEnrich requires: [Bowtie 1](http://bowtie-bio.sourceforge.net/index.shtml),
[bedtools](http://bedtools.readthedocs.org/en/latest/), 
and [samtools](http://www.htslib.org/).
RepEnrich also requires a bowtie1 indexed genome in fasta format
available. (Example `mm9.fa`) 
The RepEnrich python scripts also use [BioPython](http://biopython.org) which
can be installed with the following command:

    pip install BioPython


### Step 1) Attain repetitive element annotation
The RepEnrich setup script will build the annotation
required by RepEnrich. The default is a repeatmasker file which can be
downloaded from [repeatmasker.org](http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html),
(for instance, find the `mm9.fa.out.gz` download
[here](http://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz).
Once you have downloaded the file you can unzip it and rename it:

    gunzip mm9.fa.out.gz
    mv mm9.fa.out mm9_repeatmasker.txt

This is what the file looks like:

	SW perc perc perc query position in query matching repeat position in repeat score div. del. ins. sequence begin end (left) repeat class/family begin end  (left) ID
	687 17.4 0.0 0.0 chr1 3000002 3000156 (194195276) C L1_Mur2 LINE/L1 (4310) 1567 1413 1
	917 21.4 11.4 4.5 chr1 3000238 3000733 (194194699) C L1_Mur2 LINE/L1 (4488) 1389 913 1
	845 23.3 7.6 11.4 chr1 3000767 3000792 (194194640) C L1_Mur2 LINE/L1 (6816) 912 887 1
	621 25.0 6.5 3.7 chr1 3001288 3001583 (194193849) C Lx9 LINE/L1 (1596) 6048 5742 3

The RepEnrich setup script will also allow
you to build the annotation required by RepEnrich for a custom set of
elements using a bed file. So if you want to examine mm9 LTR repetitive
elements; you can build this file using the the repeatmasker track from
[UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables).

To do this, select genome `mm9`, click the edit box next to _Filter_, fill
in the repclass does match with `LTR`, then click submit. Back at the table
browser select option `Selected fields from primary and related tables`, 
name the output file something like `mm9_LTR_repeatmasker.bed`, and click
`Get output`. On the next page select `genoName`, `genoStart`, `genoEnd`,
`repName`, `repClass`, `repFamily` then download the file.

The UCSC puts a header on the file that needs to be removed: 

    tail -n +3 mm9_LTR_repeatmasker.bed | head -n -4 > mm9_LTR_repeatmasker_fix.bed
    mv mm9_LTR_repeatmasker_fix.bed mm9_LTR_repeatmasker.bed

This is what our custom mm9 LTR retrotransposon bed file looks like:

	$ head mm9_LTR_repeatmasker.bed

	chr1 3001722 3002005 RLTR25A LTR ERVK
	chr1 3002051 3002615 RLTR25A LTR ERVK
	chr1 3016886 3017193 RLTRETN_Mm LTR ERVK
	chr1 3018338 3018653 RLTR14 LTR ERV1

Note: It is important to get the column format right: 

* Column 1: Chromosome
* Column 2: Start
* Column 3: End
* Column 4: Class
* Column 5: Family

The file should be tab delimited. If there is no information on class
or family, you can replace these columns with the repeat name or an
arbitrary label such as `group1`.

### Step 2) Run the setup for RepEnrich

Now that we have our annotation files we can move
on to running the setup for RepEnrich. First load the dependencies
(if you use Environment Modules - otherwise just make sure that these
programs are available in your `PATH`).

    module load bowtie
    module load bedtools
    module load samtools

Next run the setup using the type of annotation you have selected (default):

    python RepEnrich_setup.py /data/mm9_repeatmasker.txt /data/mm9.fa /data/setup_folder_mm9

custom bed file:

    python RepEnrich_setup.py /data/mm9_LTR_repeatmasker.bed /data/mm9.fa /data/setup_folder_mm9 --is_bed TRUE

The previous commands have setup RepEnrich annotation that is used in
downstream analysis of data. You only have to do the setup step once for
an organism of interest. One cautionary note is that RepEnrich is only
as reliable as the genome annotation of repetitive elements for your
organism of interest. Therefore, RepEnrich performance may not be
optimal for poorly annotated genomes. 

### Step 3) Map the data to the genome using bowtie1

After the setup of the RepEnrich we now have to
map our data uniquely to the genome before running RepEnrich. This is
because RepEnrich treats unique mapping and multi-mapping reads
separately. This requires use of specific bowtie options. The bowtie
command below is recommended for RepEnrich:

    bowtie /data/mm9 -p 16 -t -m 1 -S --max /data/sampleA_multimap.fastq sample_A.fastq /data/sampleA_unique.sam

An explanation of bowtie options:

* `bowtie <bowtie_index>`
* `-p 16` - 16 cpus
* `-t` - print time
* `-m 1` - only allow unique mapping
* `-S` - output SAM
* `--max multimapping.fastq` - output multimapping reads to `multimapping.fastq`
* `unique_mapping.sam` - uniquely mapping reads

For paired-end reads the bowtie command is:

    bowtie /data/mm9 -p 16 -t -m 1 -S --max /data/sampleA_multimap.fastq -1 sample_A_1.fastq -2 sample_A_2.fastq /data/sampleA_unique.sam

The Sam file should be converted to a bam file with samtools:

    samtools view -bS sampleA_unique.sam > sampleA_unique.bam
    samtools sort sampleA_unique.bam sampleA_unique_sorted
    mv sampleA_unique_sorted.bam sampleA_unique.bam
    samtools index sampleA_unique.bam
    rm sampleA_unique.sam

You should now compute the total mapping reads for your alignment. This
includes the reads that mapped uniquely (`sampleA_unique.bam`)
and more than once (`sample_A_multimap.fastq`). The `.out` file
from your bowtie batch script contains this information (or `stdout` 
from an interactive job). 

It should looks like this:

	Seeded quality full-index
	search: 00:32:26
		# reads processed: 92084909
		# reads with at least one reported alignment: 48299773 (52.45%)
		# reads that failed to align: 17061693 (18.53%)
		# reads with alignments suppressed due to -m: 26723443 (29.02%)
	Reported 48299773 alignments to 1 output stream(s)

The total mapping reads is the `# of reads processed` - 
`# reads that failed to align`. Here our total mapping reads are:
`92084909 - 17061693 = 75023216` 


### Step 4) Run RepEnrich on the data

Now we have all the information we need to run RepEnrich.
Here is an example (for default annotation):

    python RepEnrich.py /data/mm9_repeatmasker.txt /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap.fastq sampleA_unique.bam --cpus 16

for custom bed file annotation:

    python RepEnrich.py /data/mm9_LTR_repeatmasker.bed /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap.fastq sampleA_unique.bam --is_bed TRUE --cpus 16

An explanation of the RepEnrich command:

	python RepEnrich.py 
		<repeat_annotation>
		<output_folder>
		<output_prefix>
		<RepEnrich_setup_folder>
		<multimapping_reads.fastq>
		<unique_mapping_reads.bam>
		(--is_bed TRUE)
		(--cpus 16)

If you have paired-end data the command is very similar. There will be two `sampleA_multimap.fastq` and `sampleA_multimap_1.fastq` and
`sampleA_multimap_2.fastq` from the bowtie step.

The command for running
RepEnrich in this case is (for default annotation):

    python RepEnrich.py /data/mm9_repeatmasker.txt /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap_1.fastq --fastqfile2 sampleA_multimap_2.fastq sampleA_unique.bam --cpus 16 --pairedend TRUE

for custom bed file annotation:

    python RepEnrich.py /data/mm9_LTR_repeatmasker.bed /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap_1.fastq --fastqfile2 sampleA_multimap_2.fastq sampleA_unique.bam --is_bed TRUE --cpus 16 --pairedend TRUE


### Step 5) Processing the output of RepEnrich

The final outputs will be
in the path `/data/sample_A`. This will include a few files. The most
important of which is the `sampleA_fraction_counts.txt` file. This is
the estimated counts for the repeats. I use this file to build a table
of counts for all my conditions (by pasting the individual
`*_fraction_counts.txt` files together for my complete experiment).

You can use the compiled counts file to do differential expression analysis
similar to what is done for genes. We use
[EdgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
or [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html)
to do the differential expression analysis. These are R packages that you can
download from [bioconductor](http://bioconductor.org/).

When running the EdgeR differential expression analysis
you can follow the examples in the EdgeR manual. I manually input the
library sizes (the total mapping reads we obtained in the tutorial).
Some of the downstream analysis, though, is left to your discretion.
There are multiple ways you can do the differential expression analysis.
I use the `GLM` method within the EdgeR packgage, although DESeq has
similar methods and EdgeR also has a more straightforward approach
called `exactTest`. Below is a sample EdgeR script used to do the
differential analysis of repeats for young, old, and very old mice. The
file `counts.csv` contains the ouput from RepEnrich that was made by
pasting the individual `*_fraction_counts.txt` files together for my
complete experiment. 

## Example Script for EdgeR differential enrichment analysis

```r
# EdgeR example

# Setup - Install and load edgeR
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library('edgeR')

# Input the count table and load the edgeR library.
data <- read.csv(file = "counts.csv")
library(edgeR)

#' Load the RepEnrich results - fraction counts
young_r1 <- read.delim('young_r1_fraction_counts.txt', header=FALSE)
young_r2 <- read.delim('young_r2_fraction_counts.txt', header=FALSE)
young_r3 <- read.delim('young_r3_fraction_counts.txt', header=FALSE)
old_r1 <- read.delim('old_r1_fraction_counts.txt', header=FALSE)
old_r2 <- read.delim('old_r2_fraction_counts.txt', header=FALSE)
old_r3 <- read.delim('old_r3_fraction_counts.txt', header=FALSE)
v_old_r1 <- read.delim('veryold_r1_fraction_counts.txt', header=FALSE)
v_old_r2 <- read.delim('veryold_r2_fraction_counts.txt', header=FALSE)
v_old_r3 <- read.delim('veryold_r3_fraction_counts.txt', header=FALSE)

#' Build a counts table
counts <- data.frame(
  row.names = young_r1[,1],
  young_r1 = young_r1[,4], young_r2 = young_r2[,4], young_r3 = young_r3[,4],
  old_r1 = old_r1[,4], old_r2 = old_r2[,4], old_r3 = old_r3[,4],
  v_old_r1 = v_old_r1[,4], v_old_r2 = v_old_r2[,4], v_old_r3 = v_old_r3[,4]
)

# Build a meta data object. I am comparing young, old, and veryold mice.
# I manually input the total mapping reads for each sample.
# The total mapping reads are calculated using the bowtie logs:
# # of reads processed - # reads that failed to align
meta <- data.frame(
	row.names=colnames(counts),
	condition=c("young","young","young","old","old","old","veryold","veryold","veryold"),
	libsize=c(24923593,28340805,21743712,16385707,26573335,28131649,34751164,37371774,28236419)
)

# Define the library size and conditions for the GLM
libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(meta$condition)

# Build a DGE object for the GLM
y <- DGEList(counts=counts, lib.size=libsize)

# Normalize the data
y <- calcNormFactors(y)
y$samples
plotMDS(y)

# Estimate the variance
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)

# Build an object to contain the normalized read abundance
logcpm <- cpm(y, log=TRUE, lib.size=libsize)
logcpm <- as.data.frame(logcpm)
colnames(logcpm) <- factor(meta$condition)

# Conduct fitting of the GLM
yfit <- glmFit(y, design)

# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

# Make the comparisons for the GLM
my.contrasts <- makeContrasts(
	veryold_old = veryold – old,
	veryold_young = veryold – young,
	old_young = old – young,
	levels = design
)

# Define the contrasts used in the comparisons
allcontrasts = c(
	"veryold_old",
	"veryold_young",
	"old_young"
)

# Conduct a for loop that will do the fitting of the GLM for each comparison
# Put the results into the results objects
for(current_contrast in allcontrasts) {
	lrt <- glmLRT(yfit, contrast=my.contrasts[,current_contrast])
	plotSmear(lrt, de.tags=rownames(y))
	title(current_contrast)
	res <- topTags(lrt,n=dim(c)[1],sort.by="none")$table
	colnames(res) <- paste(colnames(res),current_contrast,sep=".")
	results <- cbind(results,res[,c(1,5)])
	logfc <- cbind(logfc,res[c(1)])
}

# Add the repeat types back into the results.
# We should still have the same order as the input data
results$class <- young_r1[,2]
results$type <- young_r1[,3]

# Sort the results table by the logFC
results <- results[with(results, order(-abs(logFC.old_young))), ]

# Save the results
write.table(results, 'results.txt', quote=FALSE, sep="\t")

# Plot Fold Changes for repeat classes and types
for(current_contrast in allcontrasts) {
  logFC <- results[, paste0("logFC.", current_contrast)]
  # Plot the repeat classes
  classes <- with(results, reorder(class, -logFC, median))
  par(mar=c(6,10,4,1))
  boxplot(logFC ~ classes, data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log(Fold Change)", main=current_contrast)
  abline(v=0)
  # Plot the repeat types
  types <- with(results, reorder(type, -logFC, median))
  boxplot(logFC ~ types, data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log(Fold Change)", main=current_contrast)
  abline(v=0)
}

```


Note that the objects `logfc` contains the differential expression for the
contrast, `logcpm` contains the normalized read abundance, and `result`
contains both the differential expression and the false discovery rate for
the experimental comparison. I recommended reading more about these in the
[EdgeR manual](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).

