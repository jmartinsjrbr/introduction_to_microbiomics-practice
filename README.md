# introduction_to_microbiomics-practice

## Practices guidelines
In this section you will find a guideline for each practice that will be offered during the course

### Metagenomics: from reads to MAGs using metaWRAP pipeline
In this hands on practice you will have the opportunity to learn all the steps involved in a standard metagenomics analysis, from pre-processing steps (read trimming, host read removal and reads taxonomic assignment) to Binning steps (Binning, binning refinement, bins reassemble, Bin quantification to address their abundance, Bin classification and functional annotation).

As an example, we will use a a public dataset from NCBI SRA (https://www.ncbi.nlm.nih.gov/sra/?term=soil+metagenome+shotgun+SRR16888420).

![image](https://user-images.githubusercontent.com/11639261/142202653-eb035609-40d9-4430-8b95-29d4d8d73c61.png)
(https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR16888420).


#### 1- Download raw data and move into *DATA_FOLDER* 
```
ANALYSIS_FOLDER=$HOME/omics_course/metaWRAP_analysis
DATA_FOLDER=$HOME/omics_course/metaWRAP_data
ASSEMBLY_FOLDER=$ANALYSIS_FOLDER/ASSEMBLY_SRR
test -d $DATA_FOLDER || mkdir $DATA_FOLDER
test -d $ANALYSIS_FOLDER || mkdir $ANALYSIS_FOLDER
test -d $ASSEMBLY_FOLDER || mkdir $ASSEMBLY_FOLDER

cd $DATA_FOLDER 
```
- Download fastq files to *$HOME/omics_course/metaWRAP_data* you have created in the steps above, by clinking in the links below:
  - https://filesender.rnp.br/download.php?token=b6fe375e-af11-4eb0-a76b-6bcd96c1f1b6&files_ids=114205 
  - https://filesender.rnp.br/download.php?token=b6fe375e-af11-4eb0-a76b-6bcd96c1f1b6&files_ids=114204
  

#### 2- Run metaWRAP-Read_qc to trim the reads and remove host contamination 
```
##metaWRAP pipeline
#eval "$(conda shell.bash hook)"
#activate metawrap-env
#conda activate metawrap-env

cd $ANALYSIS_FOLDER
#mkdir 1_READ_QC
#metawrap read_qc -1 $DATA_FOLDER/ERR011347_1.fastq -2 $DATA_FOLDER/ERR011347_2.fastq -t 8 -o READ_QC/ERR011347
#metawrap read_qc -1 $DATA_FOLDER/ERR011348_1.fastq -2 $DATA_FOLDER/ERR011348_2.fastq -t 8 -o READ_QC/ERR011348
#metawrap read_qc -1 $DATA_FOLDER/ERR011349_1.fastq -2 $DATA_FOLDER/ERR011349_2.fastq -t 8 -o READ_QC/ERR011349
```
 - Pre-processing (metawrap example data)
![image](https://user-images.githubusercontent.com/11639261/142344583-556e15e7-57ca-458e-8710-f3e073d45a85.png)

 - Post-processing (metawrap example data)
 ![image](https://user-images.githubusercontent.com/11639261/142344704-5595768b-6c0a-42c0-a92b-f972b3073a6a.png)


```
##Alternatively, we will run fastp for practice purposes
CLEAN_READS=$ANALYSIS_FOLDER/1_FASTP
test -d $CLEAN_READS || mkdir $CLEAN_READS

fastp -p --thread 8 --qualified_quality_phred 15 --length_required 40 --detect_adapter_for_pe --reads_to_process 0 --low_complexity_filter --complexity_threshold 15 --report_title $CLEAN_READS/srr_report --json $CLEAN_READS/srr_fastp.json --html $CLEAN_READS/srr_fastp.html \
--in1 $DATA_FOLDER/srr16888420_1.fastq --in2 $DATA_FOLDER/srr16888420_2.fastq --out1 $CLEAN_READS/srr16888420_1.CLEAN.fastq \
--out2 $CLEAN_READS/srr16888420_2.CLEAN.fastq --unpaired1 $CLEAN_READS/srr16888420_1.CLEAN_unpaired.fastq \
--unpaired2 $CLEAN_READS/srr16888420_2.CLEAN_unpaired.fastq

##Quality repot by multiqc
eval "$(conda shell.bash hook)"
conda activate omics_course
cd $CLEAN_READS
multiqc . -o ./multiqc_report

```
 - Multiqc plot made from fasp output (metawrap example data)
 ![fastp_filtered_reads_plot](https://user-images.githubusercontent.com/11639261/142344938-294bb639-03b5-4f9a-a459-890c8409c8b7.png)

   - SRA example data (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR16888420)
   ![fastp_filtered_reads_plot (1)](https://user-images.githubusercontent.com/11639261/142345320-0fdafcd5-72dc-42d4-b978-035517093d4d.png)


#### 3- Run metaWRAP-ASSEMBLY with metaspades 
**WILL NOT BE RUN DURING DEMONSTRATION DUE TO HIGH COMPUTATIONAL COST FOR VM**
```
cat CLEAN_READS/srr*_1.CLEAN.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/srr*_2.CLEAN.fastq > CLEAN_READS/ALL_READS_2.fastq

#metawrap assembly -1 $CLEAN_READS/ALL_READS_1.fastq -2 $CLEAN_READS/ALL_READS_2.fastq -m 300 -t 16 --metaspades -o ASSEMBLY_srr
```
- **Download the final assembly file** 
```
cd $ASSEMBLY_FOLDER
wget https://github.com/jmartinsjrbr/introduction_to_microbiomics-practice/blob/main/final_assembly.fasta.gz
gunzip *.gz

```
- Assembly metrics (metawrap example data)
![Slide2](https://user-images.githubusercontent.com/11639261/142346482-924ecde4-6565-45e5-8229-333783e470cf.jpg)

- SRA example data (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR16888420)
![Slide1](https://user-images.githubusercontent.com/11639261/142346791-34c41927-b3d6-4038-b561-b8e364909ae5.jpg)

#### 4- Run metaWRAP-BINNING combining metabat2, maxbin2 and concoct 
```
conda deactivate

#activate metawrap-env
conda activate metawrap-env

###Binning###
metawrap binning -o INITIAL_BINNING_srr -t 8 -a $ASSEMBLY_FOLDER/final_assembly.fasta \
--metabat2 --maxbin2 --concoct $CLEAN_READS/ALL_READS_1.fastq $CLEAN_READS/ALL_READS_2.fastq
```

#### 5- Run metaWRAP-BINNING_REFINEMENT 
**WILL NOT BE RUN DURING DEMONSTRATION DUE TO HIGH COMPUTATIONAL COST FOR VM**

```
#metawrap bin_refinement -o BIN_REFINEMENT_srr -t 8 -A INITIAL_BINNING_srr/metabat2_bins/ -B INITIAL_BINNING_srr/maxbin2_bins/ -C INITIAL_BINNING_srr/concoct_bins/ -c 50 -x 15
```

#### 6- Run metaWRAP-REASSEMBLE_BINS
**WILL NOT BE RUN DURING DEMONSTRATION DUE TO HIGH COMPUTATIONAL COST FOR VM**

```
#metawrap reassemble_bins -o BIN_REASSEMBLY_srr -1 $CLEAN_READS/ALL_READS_1.fastq -2 $CLEAN_READS/ALL_READS_2.fastq -t 8 \
-m 300 -c 50 -x 15 -b BIN_REFINEMENT_srr/metawrap_50_15_bins

```

- **Download final bins**
``` 
#download final bins
cd $ANALYSIS_FOLDER
mkdir bins_consolidated
cd bins_consolidated

wget https://github.com/jmartinsjrbr/introduction_to_microbiomics-practice/blob/main/bin.1.permissive.fa.gz
gunzip *.gz

get https://github.com/jmartinsjrbr/introduction_to_microbiomics-practice/blob/main/binning_demo.tar.gz
tar -vzxf *.gz
rm *.gz
cd ..
```
#### 7- Run GTDBT to assign taxonomy to bins/MAGs
```
#gtdbtk classify_wf --genome_dir bins_consolidated --out_dir gtdbtk_out --prefix omics_course.gtdbtk --cpus 10 --pplacer_cpus 1 --extension fa
```

#### 8- Run PROKKA to assing annotation to MAGs
```
conda deactivate
conda activate prokka
for f in $(ls bins_consolidated/*.fa) 
do
  prokka --quiet --cpus 8 --outdir $(basename $f ".fa")_prokka_out --prefix $(basename $INFILE ".fa") --addgenes --addmrna --kingdom 'Bacteria' $INFILE

done;
 
```

#### 9- Run dbcan to CAZyme prediction
```
conda deactivate
conda activate run_dbcan
 
DB_DIR=/home/bioinformatica/omics_course/db

for f in $(ls *_prokka_out/*.faa);do
         run_dbcan.py $f protein --out_dir $(basename $f ".faa")_out --hmm_cpu 8 --tools hmmer --db_dir $DB_DIR --out_pre $(basename $f ".faa")_ --gram all
 
done

```
