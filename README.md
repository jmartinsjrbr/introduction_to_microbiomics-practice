# introduction_to_microbiomics-practice

## Practices guidelines
In this section you will find a guideline for each practice that will be offered during the course

### Metagenomics: from reads to MAGs using metaWRAP pipeline
In this hands on practice you will have the opportunity to learn all the steps involved in a standard metagenomics analysis, from pre-processing steps (read trimming, host read removal and reads taxonomic assignment) to Binning steps (Binning, binning refinement, bins reassemble, Bin quantification to address their abundance, Bin classification and functional annotation).

As an example, we will use a reduced version of the metaHIT survey dataset available at .... with metawrap pipeline. (https://github.com/bxlab/metaWRAP)

![image](https://user-images.githubusercontent.com/11639261/141855861-4383c93f-40a0-4d66-bdaa-791b33fddaaf.png)
(https://www.gutmicrobiotaforhealth.com/metahit/).


To run the analysis you can simply run the following scripts inside the VM environmet you previously prepared as described in sections above.
#### 1- Download raw data and move into *DATA_FOLDER* 
```
ANALYSIS_FOLDER=$HOME/omics_course/metaWRAP_analysis
DATA_FOLDER=$HOME/omics_course/metaWRAP_data
test -d $DATA_FOLDER || mkdir $DATA_FOLDER
test -d $ANALYSIS_FOLDER || mkdir $ANALYSIS_FOLDER
cd $DATA_FOLDER
wget https://www.dropbox.com/sh/kza1vomc5na5lo7/AADgdtGEnMww18PZMOHP-Pila?dl=0

```

#### 2- Run metaWRAP-Read_qc to trim the reads and remove human contamination 
```
##metaWRAP pipeline
#activate metawrap-env
#conda activate metawrap-env

cd $ANALYSIS_FOLDER
#mkdir 1_READ_QC
#metawrap read_qc -1 $DATA_FOLDER/ERR011347_1.fastq -2 $DATA_FOLDER/ERR011347_2.fastq -t 8 -o READ_QC/ERR011347
#metawrap read_qc -1 $DATA_FOLDER/ERR011348_1.fastq -2 $DATA_FOLDER/ERR011348_2.fastq -t 8 -o READ_QC/ERR011348
#metawrap read_qc -1 $DATA_FOLDER/ERR011349_1.fastq -2 $DATA_FOLDER/ERR011349_2.fastq -t 8 -o READ_QC/ERR011349


##Alternatively, we will run fastp for practice purposes
CLEAN_READS=$ANALYSIS_FOLDER/1_FASTP
test -d $CLEAN_READS || mkdir $CLEAN_READS

fastp -p --thread 8 --qualified_quality_phred 15 --length_required 40 --detect_adapter_for_pe --reads_to_process 0 --low_complexity_filter --complexity_threshold 15 --report_title $CLEAN_READS/srr_report --json $CLEAN_READS/srr_fastp.json --html $CLEAN_READS/srr_fastp.html \
--in1 $DATA_FOLDER/srr16888420_1.fastq --in2 $DATA_FOLDER/srr16888420_2.fastq --out1 $CLEAN_READS/srr16888420_1.CLEAN.fastq \
--out2 $CLEAN_READS/srr16888420_2.CLEAN.fastq --unpaired1 $CLEAN_READS/srr16888420_1.CLEAN_unpaired.fastq \
--unpaired2 $CLEAN_READS/srr16888420_2.CLEAN_unpaired.fastq

```

#### 3- Run metaWRAP-ASSEMBLY with metaspades 
** WILL NOT BE RUN DURING DEMONSTRATION DUE TO HIGH COMPUTATIONAL COST FOR VM **
```
cat CLEAN_READS/srr*_1.CLEAN.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/srr*_2.CLEAN.fastq > CLEAN_READS/ALL_READS_2.fastq

metawrap assembly -1 $CLEAN_READS/ALL_READS_1.fastq -2 $CLEAN_READS/ALL_READS_2.fastq -m 300 -t 16 --metaspades -o ASSEMBLY_srr
```

#### 4- Run metaWRAP-BINNING with metaspades 

### Metatranscriptomics: Metatranscriptome analysis using Sequence Annotation (SAMSA2) Pipeline 

![image](https://user-images.githubusercontent.com/11639261/141863319-65a37b17-11a6-4573-a229-d24b0d511537.png)

#### Example Files and Workflow 
For this practical demonstration we are going to use sample files provided with SAMSA 2. These files can be found in the folder: *~/omics_course/progs/samsa2/sample_files_paired-end/1_starting_files*

 - Create and enter a folder named metatranscriptomics 
```
cd omics_course 
mkdir metatranscriptomics 
cd metatranscriptomics 
```
- Copy the sample files to the metatranscriptomics folder  
```
cp -r ~/omics_course/progs/samsa2/sample_files_paired-end/1_starting_files/*.fastq . 
```
#### Setup variables pathways

- Setup pathways
```
#!/bin/bash
Setup pathways 

# VARIABLES - Set pathway for starting_location to location of samsa2  
#0. Set starting location: 
starting_location=~/omics_course/progs/samsa2 

#00. Starting files location 
starting_files_location=~/omics_course/metatranscriptomics 

#1. PEAR 
pear_location=$starting_location/programs/pear-0.9.10-linux-x86_64/bin 

#2. Trimmomatic 
trimmomatic_location=$starting_location/programs/Trimmomatic-0.36 

#3. SortMeRNA 
sortmerna_location=$starting_location/programs/sortmerna-2.1 

#4. DIAMOND 
#diamond_database="$starting_location/full_databases/RefSeq_bac" 
#diamond_subsys_db="$starting_location/full_databases/subsys_db" 
diamond_database="$starting_location/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB" 
diamond_subsys_db="$starting_location/setup_and_test/tiny_databases/subsys_db_TINY_24MB" 
diamond_location="$starting_location/programs/diamond" 

#5. Aggregation 
python_programs=$starting_location/python_scripts 
#RefSeq_db="$starting_location/full_databases/RefSeq_bac.fa" 
#Subsys_db="$starting_location/full_databases/subsys_db.fa" 
RefSeq_db="$starting_location/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB.fa" 
Subsys_db="$starting_location/setup_and_test/tiny_databases/subsys_db_TINY_24MB.fa" 

#6. R scripts and paths 
export R_LIBS="$starting_location/R_scripts/packages" 
R_programs=$starting_location/R_scripts 
```
#### Preprocessing 
The first steps in analyzing metatranscriptomics data will be the preprocessing of reads or reads QC to remove low-quality data. In SANSA pipeline, if performed paired-end sequencing, reads must be merged and filtered to remove low-quality reads and adaptor contamination. 

- Paired-end read merging 
If paired-end sequencing was performed, there will be two FASTQ sequence files for each sample provided; Sample names have “R1” and “R2”.  R1 contains the forward reads, while R2 contains the reverse reads.   

One way to merge these two files is using a read merging program, such as PEAR (Paired End reAd mergeR)( https://sco.h-its.org/exelixis/web/software/pear/ ).  
To use PEAR for merging two paired-end reads, run it from the command line, using the following command: 
```
./pear –f forward_reads.fastq –r reverse_reads.fastq  
```
 To see PEAR options, use the following command:  
```
./pear –help | less 
```
- Make a script to run all files 
```
cd $starting_files_location 
for file in $starting_files_location/*R1* 
do 
     file1=$file 
     file2=`echo $file1 | awk -F "R1" '{print $1 "R2" $2}'` 
     out_name=`echo $file | awk -F "R1" '{print $1 "merged"}'` 
     #out_name=`echo ${out_path##*/}` 
     $pear_location/pear -f $file1 -r $file2 -o $out_name 
done 
```
- To organize files create a folder and copy the merged files to the new fold 
```
mkdir $starting_files_location/step_1_output/ 
mv $starting_files_location/*merged* $starting_files_location/step_1_output/ 
```
**Result:** A fastq sequence file with overlapping paired-end reads merged for each sample.  Additionally, separate files are produced for the not combined reads; these may be included as well in the analysis. 

- Removal of adaptor contamination and low-quality reads

Raw sequences may contain low-quality reads, or reads with “contamination” – the adaptor sequences used for the process of sequencing may have been accidentally read as part of the read.  These contaminated and low-quality sequences should be removed to avoid skewing the results of a metatranscriptome analysis. 
There are many programs available for cleaning and removing adaptor sequences from raw sequence files; this pipeline is set up to use Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic).  This Java application includes primer sequences for Illumina machines, and can be used to filter out adaptor sequences from either single-end or paired-end sequencing. 
The command is structured like so: 
```
java -jar trimmomatic-0.33.jar SE –phred33 $infile $outfile_name SLIDINGWINDOW:4:15 MINLEN:99 
```
Details on these parameters, as well as other commands, can be found in the Trimmomatic manual.  That manual can be accessed here: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf   
**Result:** A fastq sequence file with low-quality sequences and adaptor contamination removed. 

```
java -jar $trimmomatic_location/trimmomatic-0.36.jar SE -phred33 $starting_files_location/step_1_output/control_1_TINY_merged.assembled.fastq $starting_files_location/step_1_output/control_1_TINY_clean.fastq SLIDINGWINDOW:4:15 MINLEN:99 
```

- Make a script to run all files 
```
for file in $starting_files_location/step_1_output/*.assembled.* 
do 
     name=`echo $file | awk -F "merged" '{print $1 "clean.fastq"}'` 
     java -jar $trimmomatic_location/trimmomatic-0.36.jar SE -phred33 $file $name SLIDINGWINDOW:4:15 MINLEN:99 
done 
```
- To organize files create a folder and copy the merged files to the new fold 
```
mkdir $starting_files_location/step_2_output/ 
mv $starting_files_location/step_1_output/*clean.fastq $starting_files_location/step_2_output/ 
```
- GETTING RAW SEQUENCES COUNTS 
```
for file in $starting_files_location/step_2_output/*clean.fastq 
do 
     python3 $python_programs/raw_read_counter.py -I $file -O $starting_files_location/step_2_output/raw_counts.txt 
done 
```

- Removal of ribosome sequences 
One issue when sequencing extracted RNA is filtering out mRNA from the much more common ribosomal RNA, or rRNA. Although rRNA comprises the majority of all extracted RNA from microbiome communities, it can obscure the more important mRNAs.  For best results, ribodepletion methods should be used on the biological samples after RNA extraction and before sequencing as a quality control step. 

However, biological ribodepletion kits are not completely effective at removing all ribosomes.  Stripping out remaining ribosomal reads will help increase the speed of downstream pipeline steps, by resulting in fewer total reads to be annotated (the slowest and most computationally intensive step) and analyzed. 

SortMeRNA (http://bioinfo.lifl.fr/RNA/sortmerna/) is a robust ribosomal read filtering tool that can incorporate multiple databases (SILVA, GreenGenes, RDP) for rRNA identification.  Note that SortMeRNA was originally designed to select rRNA sequences, rather than to remove them, and so the reads discarded by SortMeRNA are, in fact, the mRNAs needed for metatranscriptome analysis. 
For detailed instructions on using SortMeRNA, be sure to consult the included user manual for version 2.1.   

Note that the “--other" flag MUST be applied when using SortMeRNA!  Without this flag, all reads that do not match the ribosomal RNAs in the reference will be discarded.  These reads are, in fact, the mRNAs, and must be preserved for the following steps in the SAMSA pipeline. 

Example SortMeRNA command (matching against the 16S SILVA bacterial database, included in SortMeRNA download): 
```
sortmerna --refsilva-bac-16s-db --reads $file.fastq --aligned $file.ribosomes --other $file.ribodepleted --fastx --num_alignments 1 --log –v 
```
From this command, two files will be produced; the $file.ribosomes will contain all sequences from the original file identified as rRNA, while the $file.ribodepleted will contain all reads discarded by SortMeRNA (aka not identified as ribosomes, to be used in the next step of the SAMSA pipeline). 

Note that while the identified ribosomal sequences can potentially be used for other analyses, the biological ribodepletion that is strongly recommended for all samples before sequencing will likely skew these results, making them unusable for organism-specific abundance measurements. 

**Result:** A fastq sequence file with ribosomal sequences removed; additionally, a second file is created containing said ribosomal sequences for optional taxonomic profiling. 
```
#REMOVING RIBOSOMAL READS WITH SORTMERNA 
# Note: this step assumes that the SortMeRNA databases are indexed. If not, do that first (see the SortMeRNA user manual for details). 
for file in $starting_files_location/step_2_output/*clean.fastq 
do 
  name=`echo $file | awk -F "clean" '{print $1 "ribodepleted"}'` 
sortmerna_location/sortmerna \ 
--ref $sortmerna_location/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/index/silva-bac-16s-db \ 
--reads $file \ 
--aligned $file.ribosomes \ 
--other $name \ 
--fastx \ 
--num_alignments 0 \ 
--log -v 
done 

mkdir $starting_files_location/step_3_output/ 

mv $starting_files_location/step_2_output/*ribodepleted* $starting_files_location/step_3_output/ 
```

#### Annotation 
Now that the initial sequences have been merged (if using paired-end sequencing), cleaned, and stripped of adaptor sequence contamination and ribosomal sequences, the next step is to annotate the mRNA reads to their corresponding match in a reference database.  

It can be done using DIAMOND, a superfast BLAST-like aligner tool. This tool can process reads up to 10,000x as fast as BLASTX, with very little loss in accuracy. DIAMOND can also annotate against any provided database, allowing for custom databases to be created and searched against. 

- Creating a DIAMOND-indexed database 
Any database file needs to be indexed by DIAMOND and converted into a binary file before it can be searched against.  DIAMOND will convert any fasta file to a usable database with the following command: 

```
diamond makedb --in $database --db $database 
```

The --in flag specifies the starting fasta file that will be converted to a DIAMOND-searchable database. 

DIAMOND can be given any database file to be indexed. Two databases that should be interesting to microbiome researchers are the NCBI RefSeq database and the SEED Subsystems database.  Maintained by NCBI, RefSeq is one of the most complete databases for general purposes and is generally accepted to contain high-quality annotations.  SEED Subsystems offers the unique ability to sort specific functions into hierarchies, letting similar functions be grouped under a category heading, such as “cellular respiration” or “protein biosynthesis.”  This can be very useful for examining overall functional activity within a metatranscriptome. 

The RefSeq database can be accessed through NCBI’s FTP site, here: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/ .  The simplest approach is to download all non-redundant protein sequence files, use cat on the command line to merge them together into a single gzipped file: 

```
cat file1.gz file2.gz file3.gz > all_files.gz 
```

The SEED Subsystems database can be accessed through their FTP site here: ftp://ftp.theseed.org/subsystems/ 

Note, however, that the SEED Subsystems database is not readily downloadable in a fasta format that can be indexed by DIAMOND.  The different levels of Subsystems hierarchy are maintained in different files.   For merging these files together to create a single, indexable database that contains all hierarchy information, see the relevant Github repository here: https://github.com/transcript/subsystems  

- Annotating a file against a DIAMOND database 
```
For demonstration purpose we will use a TINY version of RefSeq and Subsystems database. 
for file in $starting_files_location/step_3_output/*ribodepleted.fastq 
do 
     name=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'` 

     $diamond_location blastx --db $diamond_database -q $file -f 6 -o $name -k 1 --sensitive 
done 

mkdir $starting_files_location/step_4_output/ 
mv $starting_files_location/step_3_output/*annotated* $starting_files_location/step_4_output/ 

diamond blastx --db $diamond_database -q $filename -f 6 -o $diamond_output  -k 1  
```
The resulting data table is ready for step 3 in the SAMSA pipeline: **aggregation**. 

#### Aggregation 

Now that each metatranscriptome file has been annotated, the next step is to reduce the results down into a condensed and simplified format for statistical analysis.  DIAMOND returns the best match for each read in the starting file that meets its parameters for sequence specificity, much like a line-item receipt from a grocery store.  This step converts this large file into a condensed, sorted summary table that returns the total number of hits to each specific organism or function. 

NOTE: This step will create two summary files for each starting metatranscriptome; one file will contain annotations grouped by organism (all Bacteroides reads will be grouped together), while the other file will contain annotations grouped by function (all reads coding for the enzyme lactase will be grouped together).  Later steps will document the steps necessary to perform a search for all functions expressed by a specific organism or group of organisms, or vice versa, all organisms performing a specific function or set of functions. 

This next step uses the Python program “DIAMOND_analysis_counter.py”, and for the RefSeq database, will require access to the original (readable, not DIAMOND-converted) database file. 

To use this program for aggregating all reads by organism: 
```
python DIAMOND_analysis_counter.py –I $infile –D database_file –O 
```
 
And to use this program for aggregating all reads by function: 
```
python DIAMOND_analysis_counter.py –I $infile –D database_file –F 
``` 
The result is a 3 column table, saved in tab-separated values (.tsv) format.  The columns are as follows: 

*(percentage of total reads)		(read count)		(annotated organism or function)*

The resulting files can either be viewed directly, or can be imported into R for further statistical analysis and figure generation. 


- AGGREGATING WITH ANALYSIS_COUNTER 
```
for file in $starting_files_location/step_4_output/*RefSeq_annotated* 
do 
  python3 $python_programs/standardized_DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O 
  python3 $python_programs/standardized_DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F 
done 

mkdir $starting_files_location/step_5_output/ 
mkdir $starting_files_location/step_5_output/RefSeq_results/ 
mkdir $starting_files_location/step_5_output/RefSeq_results/org_results/ 
mkdir $starting_files_location/step_5_output/RefSeq_results/func_results/ 

mv $starting_files_location/step_4_output/*organism.tsv $starting_files_location/step_5_output/RefSeq_results/org_results/ 
mv $starting_files_location/step_4_output/*function.tsv $starting_files_location/step_5_output/RefSeq_results/func_results/ 

```

#### Analysis 
The aggregated results files generated by the Python script in step 3, aggregation, can be used for several different comparisons and analyses in R.  Note that, while some of these scripts can be run in base (command line) R, they can also be handled using RStudio, which allows for evaluation of intermediate data tables and the ability to better isolate any issues that may arise due to differences in database and output structure. 

Before using the Rscripts provided with samsa2 pipeline it is necessary to install some R packages 

```
R -e 'Install.packages(c “pheatmap”, “vegan”, “plyr”)'

# Need installation of EnhancedVolcano package beforehand 
R -e 'if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('EnhancedVolcano')'
#quit() 
```

- Statistical analysis and differential expression using DESeq2 

The R script “run_DESeq_stats.R” evaluates which organisms or functions exhibit a significantly different level when comparing between control and experimental results? 

This R script, if run from the command line, takes two arguments: the directory location of the summary files that will be analyzed as ARGV1, and the name under which the results should be saved as ARGV2.  In the specified directory, control file names should begin with “control_”, and experimental file names should begin with “experimental_”. 

**Result:** This R program outputs a tab-delimited results table that, for each organism or function in the infiles, prints out the controlMean, the experimentalMean, the log2FoldChange, and the absolute and multiple-hypothesis-adjusted p-values.   

- R ANALYSIS 
```
Rscript $R_programs/run_DESeq_stats.R -I $starting_files_location/step_5_output/RefSeq_results/org_results/ -O RefSeq_org_DESeq_results.tab -R $starting_files_location/step_2_output/raw_counts.txt 

Rscript $R_programs/run_DESeq_stats.R -I $starting_files_location/step_5_output/RefSeq_results/func_results/ -O RefSeq_func_DESeq_results.tab -R $starting_files_location/step_2_output/raw_counts.txt
```

- Diversity measures using organism data 

Two statistics, Shannon and Simpson diversity, can provide a useful overview of the diversity of a microbiome.  Two R programs evaluate this measure: “diversity_stats.R” will calculate the average diversity seen in each sample group, while “diversity_graphs.R” will create graphs of the diversity statistics for each individual sample. 

```
Rscript $R_programs/diversity_stats.R -I $starting_files_location/step_5_output/RefSeq_results/org_results/ 
Rscript $R_programs/diversity_graphs.R -I $starting_files_location/step_5_output/RefSeq_results/org_results/ 

evince $starting_files_location/step_5_output/RefSeq_results/org_results/ diversity_graph.pdf 
```

**Result:** “diversity_stats.R” will return the mean Shannon and Simpson diversities for both the control and experimental groups.   

*diversity_graphs.R* will create two graphs, one showing Shannon diversity calculated for each individual data sample, the other showing Simpson diversity calculated for each individual data sample. 

 ![image](https://user-images.githubusercontent.com/11639261/141871991-87ce8b43-f9fe-4317-85a6-f861c7f1dea4.png)

- Combined graphs of top (most abundant) organisms or annotations 

A stacked graph showing the top organisms (or functions) within all metatranscriptomes in a project provides a useful overview of results. “make_combined_graphs.R” is an R script that creates these graphs for either the top most abundant organisms or functions.  Similar to the above R programs, it reads in summaries of all metatranscriptome files and generates a stacked bar graph from this data. 

Rscript $R_programs/make_combined_graphs.R -D $starting_files_location/step_5_output/RefSeq_results/org_results/ 
evince $starting_files_location/step_5_output/RefSeq_results/org_results/ combined_graph.pdf 

**Result:** “make_combined_graphs.R” will generate a PDF containing two stacked bar graphs; the top bar graph will reveal relative expression (in percentages), while the lower bar graph will reveal absolute expression (count numbers).  A legend is displayed beneath.  Control samples are displayed on the left, while experimental samples are displayed on the right.  Each bar represents a single metatranscriptome input file. 

![image](https://user-images.githubusercontent.com/11639261/141872141-dcbcf631-7f35-493b-a5bb-e96d9218dc4b.png)

- PCA plots of metatranscriptomes  

A PCA plot is useful for showing clustering between different metatranscriptomes, demonstrating that different samples can be regarded as different groups or categories.  The R program “make_DESeq_PCA.R” uses the DESeq2 package to generate the data structure for this graph, and then creates it with ggplot2. 

```
Rscript $R_programs/make_DESeq_PCA.R -I $starting_files_location/step_5_output/RefSeq_results/org_results/ 

evince $starting_files_location/step_5_output/RefSeq_results/org_results/PCA_plot.tab.pdf 
```
 
**Result:**  A PCA plot showing experimental vs. control samples. Further customization can be performed in the script if desired to provide different shapes for different individual samples. 

![image](https://user-images.githubusercontent.com/11639261/141872251-eea0a1cf-77d5-4b73-9203-533615fb2e67.png)

#### Analysis with Subsystems database 
- Subsystems Hierarchy Comparisons 

While the detailed information of a metatranscriptome is incredibly powerful, it can feel overwhelming to consider hundreds of thousands of different annotated functions, when multiple thousands may be noted as significantly differentially expressed.  One unique method for simplifying the “at-a-glance” metatranscriptome is to use a hierarchical database, such as TheSeed’s Subsystems system of annotations. 

TheSeed’s Subsystems database (referred to hereafter as “Subsystems”) features four increasing levels of hierarchy for each annotated function.  They’re best imagined as a fractal tree, with specific functions forming the leaves.  For example, a specific sequence might be linked to TcuB tricarballylate oxidation (level 4 subsystem).  This, in turn, falls under the category of “Tricarballylate utilization” (level 3 subsystem).  This falls into the level 2 subsystem of “organic acids”, which is in the level 1 subsystem category of “Carbohydrates.” 

Another way of viewing it: 

1- Carbohydrates 

  a. Organic acids 

   i. Tricarballylate utilization 

   1. TcuB tricarballylate oxidation to cis-aconitate 

The Subsystems database can be reworked into a SAMSA-searchable database, and metatranscriptome hits against this database can be compared at each hierarchy level.  Thus, instead of sorting through 100,000+ specific functions, it’s possible to obtain a metatranscriptome “overview” by looking only at differences in the top 20 or so level 1 Subsystems categories. 

- Annotating against a flattened Subsystems database 

Unfortunately, the Subsystems database, although available for download through an FTP portal (ftp://ftp.theseed.org/subsystems/), is not pre-configured in a format that makes it indexable by algorithms like DIAMOND.  For a more detailed view of obtaining a flattened, DIAMOND-indexable version of the Subsystems database, see the following Github page: https://github.com/transcript/subsystems.  Using the Python program “subsys_db_rebuilder.py” and the other files mentioned in the Github readme, the Subsystems database can be flattened into a single file containing all hierarchy, which can then be indexed by DIAMOND and used as a searchable database, much like the NCBI RefSeq database used in earlier steps. 

Once the Subsystems database is converted to a flattened format, running a DIAMOND annotation search against it is nearly identical to searching against RefSeq, in the earlier section of this guide: 

- ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS 
```
for file in $starting_files_location/step_3_output/*ribodepleted.fastq 
do 
   name=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'` 
   $diamond_location blastx --db $diamond_subsys_db -q $file -o $name -k 1 --sensitive 
done 

mv $starting_files_location/step_3_output/*subsys_annotated* $starting_files_location/step_4_output/ 
```

- Converting annotations to hierarchies 

The DIAMOND results will return the match for each given read to the specific Fig ID, the base unit of reference in the Subsystems database.  To obtain the full hierarchy, a Python program needs to go back to the database and recheck each read against the reference, adding in the necessary information. 

The Python program that accomplishes this is called “DIAMOND_subsystems_analysis_counter.py”, and operates very similarly to the “DIAMOND_analysis_counter.py” program used for RefSeq results.  Similarly, this analysis counter program reads back through the DIAMOND output and condenses the reads down into a sorted abundance list – but it also includes all hierarchy information, if available, printed in additional columns in the tab-delimited output. 

**Usage:**
   *python DIAMOND_subsystems_analysis_counter.py -I $file -D $subsystems_db -O $file.analyzed -P $file.receipt*

Note that, as the Python program does not read binary files, the linked Subsystems database must be the flattened file created for conversion to a DIAMOND-indexable database (but NOT the converted .daa file!).  The output will be saved under the name given with the –O flag, while the optional –P flag will create another file containing a line-by-line “receipt” with all hierarchy included. 


- PYTHON SUBSYSTEMS ANALYSIS COUNTER 

```
for file in $starting_files_location/step_4_output/*subsys_annotated 
do 
  python3 $python_programs/DIAMOND_subsystems_analysis_counter.py -I $file -D $Subsys_db -O $file.hierarchy -P $file.receipt 
  python3 $python_programs/subsys_reducer.py -I $file.hierarchy 
done 

mkdir $starting_files_location/step_5_output/Subsystems_results/ 
mkdir $starting_files_location/step_5_output/Subsystems_results/receipts/ 

mv $starting_files_location/step_4_output/*.reduced $starting_files_location/step_5_output/Subsystems_results/ 
mv $starting_files_location/step_4_output/*.receipt $starting_files_location/step_5_output/Subsystems_results/receipts/ 

rm $starting_files_location/step_4_output/*.hierarchy 
```

The resulting output file from this is ready to be imported into R for analysis, either for calculating significantly differentially expressed categories of functions, or for the creation of pie charts reflecting the relative levels of reads towards each Subsystems category. 

- Analyzing Subsystems annotations for differential expression in R 
Because the Subsystems results from the Python aggregating program contain additional columns with hierarchy information, a separate R script is needed to analyze the results.  These two R scripts are named “Subsystems_DESeq_stats.R” and “Subsystems_pie_charts.R”.   

For figuring out significantly differentially expressed genes at each of the different hierarchy levels, “Subsystems_DESeq_stats.R” loads in the Subsystems results files, parses the resulting data tables down to only the requested hierarchy level, and then runs DESeq2 to figure out differential expression.   
```
Rscript $R_programs/Subsystems_DESeq_stats.R -I $starting_files_location/step_5_output/Subsystems_results/ -O Subsystems_level-1_DESeq_results.tab -L 1 -R $starting_files_location/step_2_output/raw_counts.txt 
```
**Result:** A data table containing DESeq evaluated differentially expressed transcripts for the selected hierarchy level, saved in a tab-delimited format. 

- Analyzing Subsystems annotations for creating pie charts in R 

Similar to the above program, “Subsystems_pie_charts.R” will load in the Subsystems files, compress down to only the selected level of hierarchy, and will then create a pie chart.  Note that currently, “Subsystems_pie_charts.R” does NOT choose to make two different pie charts for experimental vs. control files.  Instead, it will create a single pie chart.  If two contrasting pie charts are desired, run the program twice, feeding it each different group of files (once with the experimental files, a second time, separately, with control files). 

By default, this script will produce output graphs of level 1 results.  For different levels, change the level indicator on line 49 (note that higher levels may provide too many files and will require the use of additional colors). 

```
Rscript $R_programs/Subsystems_pie_charts.R -I $starting_files_location/step_5_output/Subsystems_results/ -L 1  
```
 
**Result:** A pie chart, saved as a PDF under the given save filename.  Details can be altered by changing the ggplot2 graph created at the end of the program. 
![image](https://user-images.githubusercontent.com/11639261/141872840-30b94bb2-df2d-40a9-910c-638c3b65a1a5.png)

#### Further Applications Functions by a specific Organism 

Although looking at the functional activity of all organisms can provide a useful overview of a metatranscriptome and may offer suggestions as to where to focus further investigation, it’s also important to have the ability to narrow the focus and only examine the outputs linked to specific organisms of interest.  For example, a gut microbiome study may be particularly interested in the activity of Lactobacillus species, and would like to restrict functional annotations to solely this genus. 

Fortunately, because each read processed by SAMSA receives both an organism and a functional annotation, this is possible!  This is made possible using a Python program named “DIAMOND_specific_organism_retreiver.py”. 

This program takes three inputs – the input file, which is the DIAMOND results file from annotating against the RefSeq database (-I flag), the specific organism to be selected (-SO flag), and the RefSeq database file (-D flag)(note: make sure to specify the original fasta file, not the DIAMOND-indexed binary).  The program first reads through the RefSeq database, constructing a dictionary of RefSeq IDs and their matching organisms.  It next reads through the infile, checking each annotated RefSeq ID’s entry to see whether the organism of interest’s name is present in the dictionary entry.  Only those entries with the chosen organism name present are passed on to the outfile (which is automatically generated, named after the combination of the target infile and the specific organism).   

**Result:** The program produces an output, named after the infile and the specific organism it’s searching against.  This output file contains all individual transcripts that originated from the specific organism; the next step in the pipeline is to use DIAMOND_analysis_counter.py to reduce this down to a sorted abundance list. 
