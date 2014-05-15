Tutorial for RepEnrich -- By Steven Criscione -- Email:
steven\_criscione@brown.edu This example is for mouse genome mm9. Before
getting started you should make sure you have installed the dependencies
for RepEnrich. RepEnrich requires bowtie1, bedtools, and samtools.
RepEnrich also requires a bowtie1 indexed genome in fasta format
available. (Example mm9.fa) 

----------

----------
**Step 1) Attain repetitive element
annotation** The RepEnrich setup script will build the annotation
required by RepEnrich. The default is a repeatmasker file which can be
downloaded from repeatmasker.org at
[http://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz](http://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz "http://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz").
Once you download the file you can unzip it and rename it.

    gunzip mm9.fa.out.gz

    mv mm9.fa.out mm9_repeatmasker.txt

This is what the file looks like:

    head mm9_repeatmasker.txt

SW perc perc perc query position in query matching repeat position in repeat score div. del. ins. sequence begin end (left) repeat class/family begin end  (left) ID<br>
687 17.4 0.0 0.0 chr1 3000002 3000156 (194195276) C L1\_Mur2 LINE/L1 (4310) 1567 1413 1<br>
917 21.4 11.4 4.5 chr1 3000238 3000733 (194194699) C L1\_Mur2 LINE/L1 (4488) 1389 913 1<br>
845 23.3 7.6 11.4 chr1 3000767 3000792 (194194640) C L1\_Mur2 LINE/L1 (6816) 912 887 1<br>
621 25.0 6.5 3.7 chr1 3001288 3001583 (194193849) C Lx9 LINE/L1 (1596) 6048 5742 3<br>

The RepEnrich setup script will also allow
you to build the annotation required by RepEnrich for a custom set of
elements using a bed file. So if you want to examine mm9 LTR repetitive
elements; you can build this file using the the repeatmasker track from
UCSC genome table browser at
[http://genome.ucsc.edu/cgi-bin/hgTables](http://genome.ucsc.edu/cgi-bin/hgTables "http://genome.ucsc.edu/cgi-bin/hgTables").
To do so, you would want to select genome mm9, click the edit box next
to filter, fill in the repclass does match with "LTR", then click
submit, back at the table browser select option "selected fields from
primary and related tables", name the output file something like
"mm9\_LTR\_repeatmasker.bed", and click get output. On the next page
select "genoName, genoStart, genoEnd, repName, repClass, repFamily" then
download the file. The UCSC puts a header on the file that can be
removed.

    tail -n +3 mm9_LTR_repeatmasker.bed | head -n -4 > mm9_LTR_repeatmasker_fix.bed

    mv mm9_LTR_repeatmasker_fix.bed mm9_LTR_repeatmasker.bed

This is what our custom mm9 LTR retrotransposon bed file looks like:

    head mm9_LTR_repeatmasker.bed

chr1 3001722 3002005 RLTR25A LTR ERVK<br> 
chr1 3002051 3002615 RLTR25A LTR ERVK<br> 
chr1 3016886 3017193 RLTRETN\_Mm LTR ERVK<br> 
chr1 3018338 3018653 RLTR14 LTR ERV1<br> 

Note: The format is
important to get accurate information on repeat name, class, and family.
Column 1: chromosome, column 2: start, column 3: end, column 4: class,
column 5: family. The file should be tab delimited. If there is no
information on class or family, you can replace these columns with the
repeat name or an arbitrary label such as "group1". 

----------

----------

**Step 2) Run the
setup for RepEnrich**<br>

Now that we have our annotation files we can move
on to running the setup for RepEnrich. First load the dependencies.

    module load bowtie

    module load bedtools

    module load samtools

Next run the setup using the type of annotation you have selected:
default:

    python RepEnrich_setup.py /data/mm9_repeatmasker.txt /data/mm9.fa /data/setup_folder_mm9

custom bed file:

    python RepEnrich_setup.py /data/mm9_LTR_repeatmasker.bed /data/mm9.fa /data/setup_folder_mm9 --is_bed TRUE

The previous commands have setup RepEnrich annotation that is used in
downstream analysis of data. You only have to do the setup step once for
an organism of interest. One cautionary note is that RepEnrich is only
as reliable as the genome annotation of repetitive elements for your
organism of interest. Therefore, RepEnrich performance may not be
optimal for poorly annotated genomes. 

----------

----------

**Step 3) Map the data to the
genome using bowtie1** 

After the setup of the RepEnrich we now have to
map our data uniquely to the genome before running RepEnrich. This is
because RepEnrich treats unique mapping and multi-mapping reads
separately. This requires use of specific bowtie options. The bowtie
command below is recommended for RepEnrich:

    bowtie /data/mm9 -p 16 -t -m 1 -S --max /data/sampleA_multimap.fastq sample_A.fastq /data/sampleA_unique.sam

An explanation of bowtie options bowtie bowtie\_index (-p 16)=16\_cpus
(-t)=print\_time (-m 1)=only\_allow\_unique\_mapping ( -S )=output\_sam
(--max)=output\_multimapping\_reads multimapping.fastq input.fastq
unique\_mapping.sam For paired-end reads the bowtie command is:

    bowtie /data/mm9 -p 16 -t -m 1 -S --max /data/sampleA_multimap.fastq -1 sample_A_1.fastq -2 sample_A_2.fastq /data/sampleA_unique.sam

The Sam file should be converted to a bam file with samtools:

    samtools view -bS sampleA_unique.sam > sampleA_unique.bam

    samtools sort sampleA_unique.bam sampleA_unique_sorted

    mv sampleA_unique_sorted.bam sampleA_unique.bam

    samtools index sampleA_unique.bam

    rm sampleA_unique.sam

You should now compute the total mapping reads for your alignment. This
includes the reads that mapped uniquely (in your sampleA\_unique.bam)
and more than once ( in your sample\_A\_multimap.fastq). The .out file
from your bowtie batch script contains this information (or stdout from
interactive bowtie running). 

It looks like:

Seeded quality full-index\n
search: 00:32:26 \# reads processed: 92084909<br>
\# reads with at least one reported alignment: 48299773 (52.45%) <br>\# reads that failed to align:
17061693 (18.53%) <br>\# reads with alignments suppressed due to -m:
26723443 (29.02%) <br>Reported 48299773 alignments to 1 output stream(s)

The total mapping reads is the (\# of reads processed) - (\# reads that
failed to align)<br>
Here our total mapping reads is: 92084909 - 17061693 =
75023216 

----------

----------


**Step 4) Run RepEnrich on the data** 

Now we have all the information we need to run RepEnrich. Here is an example: for default
annotation:

    python RepEnrich.py /data/mm9_repeatmasker.txt /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap.fastq sampleA_unique.bam --cpus 16

for custom bed file annotation:

    python RepEnrich.py /data/mm9_LTR_repeatmasker.bed /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap.fastq sampleA_unique.bam --is_bed TRUE --cpus 16

An explanation of RepEnrich command: python RepEnrich.py
repeat\_annotation output\_folder output\_prefix
RepEnrich\_setup\_folder multimapping\_reads.fastq
unique\_mapping\_reads.bam total\_mapping\_reads (--is\_bed
TRUE)=repeat\_annotation\_is\_bedfile ( --cpus 16)=use\_16\_cpus

If you have paired-end data the command is very similar. There will be two
sampleA\_multimap.fastq files sampleA\_multimap\_1.fastq and
sampleA\_multimap\_2.fastq from the bowtie step. The command for running
RepEnrich in this case is: for default annotation:

    python RepEnrich.py /data/mm9_repeatmasker.txt /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap_1.fastq --fastqfile2 sampleA_multimap_2.fastq sampleA_unique.bam 75023216 --cpus 16 --pairedend TRUE

for custom bed file annotation:

    python RepEnrich.py /data/mm9_LTR_repeatmasker.bed /data/sample_A sample_A /data/hg19_setup_folder sampleA_multimap_1.fastq --fastqfile2 sampleA_multimap_2.fastq sampleA_unique.bam 75023216 --is_bed TRUE --cpus 16 --pairedend TRUE


----------

----------

**Step 5) Processing the output of RepEnrich** 

The final outputs will be
in the path /data/sample\_A. This will include a few files. The most
important of which is the sampleA\_fraction\_counts.txt file. This is
the estimated counts for the repeats. I use this file to build a table
of counts for all my conditions (by pasting the individual
\*\_fraction\_counts.txt files together for my complete experiment). You
can use the compiled counts file to do differential expression analysis
similar to what is done for genes. We use EdgeR or DESEQ to do the
differential expression analysis. These are R packages that you can
download from bioconductor.

 When running the EdgeR differential
expression analysis
([http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html](http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html "http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html"))
you can follow the examples in the EdgeR manual. I manually input the
library sizes (the total mapping reads we obtained in the tutorial).
Some of the downstream analysis, though, is left to your discretion.
There are multiple ways you can do the differential expression analysis.
I use the GLM method within the EdgeR packgage, although DESeq has
similar methods and EdgeR also has a more straightforward approach
called exactTest. Below is a sample EdgeR script used to do the
differential analysis of repeats for young, old, and very old mice. The
file counts.csv contains the ouput from RepEnrich that was made by
pasting the individual \*\_fraction\_counts.txt files together for my
complete experiment. 

## Example Script for EdgeR differential enrichment analysis<br>
    # EdgeR example<br>
    # Input the count table and load the edgeR library.<br>
    data <- read.csv(file = “counts.csv”)<br>
    library(edgeR)<br>
    # Define counts and groups. My counts table has repeat name in the first column so I remove that from counts object and make it the rownames.<br>
    counts <- data[, -1 ]<br>
    rownames(counts) <- data[,1]<br>
    # Build a meta data object. I am comparing young, old, and veryold mice. I manually input the total mapping reads for each sample.<br>
    meta <- data.frame(<br>
    row.names=colnames(counts),<br>
    condition=c(“young”,”young”,”young”,”old”,”old”,”old”,”veryold”,”veryold”,”veryold”),<br>
    libsize=c(24923593,28340805,21743712,16385707,26573335,28131649,34751164,37371774,28236419))<br>
    # Define the library size and conditions for the GLM<br>
    libsize <- meta$libsize<br>
    condition <- factor(meta$condition)<br>
    design <- model.matrix(~0+condition)<br>
    colnames(design) <- levels(meta$condition)<br>
    # Build a DGE object for the GLM<br>
    y <- DGEList(counts=counts, lib.size=libsize)<br>
    # Normalize the data<br>
    y <- calcNormFactors(y)<br>
    # Estimate the variance<br>
    y <- estimateGLMCommonDisp(y, design)<br>
    y <- estimateGLMTrendedDisp(y, design)<br>
    y <- estimateGLMTagwiseDisp(y, design)<br>
    # Build an object to contain the normalized read abundance<br>
    logcpm = cpm(y, log=TRUE, lib.size=libsize)<br>
    logcpm = as.data.frame(logcpm)<br>
    colnames(logcpm) = factor(meta$condition)<br>
    # Conduct fitting of the GLM<br>
    yfit <- glmFit(y, design)<br>
    # Initialize result matrices to contain the results of the GLM<br>
    results = matrix(nrow=dim(counts)[1],ncol=0)<br>
    logfc = matrix(nrow=dim(counts)[1],ncol=0)<br>
    # Make the comparisons for the GLM<br>
    my.contrasts <- makeContrasts(<br>
    veryold_old = veryold – old,<br>
    veryold_young = veryold – young,<br>
    old_young = old – young ,<br>
    levels=design)<br>
    # Define the contrasts used in the comparisons<br>
    allcontrasts = c(<br>
    “veryold_old” ,<br>
    “veryold_young” ,<br>
    “old_young”<br>
    )<br>
    # Conduct a for loop that will do the fitting of the GLM for each comparison<br>
    # Put the results into the results objects<br>
    for(current_contrast in allcontrasts) {<br>
    lrt <- glmLRT(yfit, contrast=my.contrasts[,current_contrast])<br>
    res <- topTags(lrt,n=dim(c)[1],sort.by=”none”)$table<br>
    colnames(res) <- paste(colnames(res),current_contrast,sep=”.”)<br>
    results <-cbind(results,res[,c(1,5)])<br>
    logfc <-cbind(logfc,res[c(1)])<br>
    }<br>


----------

Note the objects logfc contains the differential expression for the contrast, logcpm contains the normalized read abundance, and result contains both the differential expression and the false discovery rate for the experimental comparison. I recommended reading more about these in the EdgeR manual.

