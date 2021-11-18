# introduction_to_microbiomics-practice

## Practices guidelines
In this section you will find a guideline for each practice that will be offered during the course

### Metagenomics: from reads to MAGs using metaWRAP pipeline
In this hands on practice you will have the opportunity to learn all the steps involved in a standard metagenomics analysis, from pre-processing steps (read trimming, host read removal and reads taxonomic assignment) to Binning steps (Binning, binning refinement, bins reassemble, Bin quantification to address their abundance, Bin classification and functional annotation).

As an example, we will use a a public dataset from NCBI SRA (https://www.ncbi.nlm.nih.gov/sra/?term=soil+metagenome+shotgun+SRR16888420).

![image](https://user-images.githubusercontent.com/11639261/142202653-eb035609-40d9-4430-8b95-29d4d8d73c61.png)
(https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR16888420).


-Note: **In addition, we will show as demonstration some results from megaHIT human gut survey retrieved from metaWRAP tutorial (https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md)**

 ![Slide1](https://user-images.githubusercontent.com/11639261/142384996-21dfc2a3-a1d9-4152-8679-8e6c20a03ef4.jpg)

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
- Download fastq files to *$HOME/omics_course/metaWRAP_data* you have created in the steps above, by cliking in the links below:
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
- Binning refinement (metawrap example data)
```
bin	completeness	contamination	GC	lineage	N50	size	binner
bin.5	99.32	1.006	0.406	Clostridiales	70006	2090427	binsA
bin.9	98.4	1.6	0.311	Euryarchaeota	12258	1688675	binsC
bin.6	85.42	1.027	0.374	Clostridiales	4985	2055629	binsC
bin.10	84.77	2.516	0.445	Clostridiales	2839	1452931	binsC
bin.12	79.48	1.034	0.464	Bacteria	12851	4215790	binsA
bin.7	75.33	0.459	0.287	Clostridiales	15087	1460311	binsA
bin.1	71.23	2.314	0.422	Selenomonadales	3450	1658218	binsC
bin.4	70.90	1.995	0.509	Clostridiales	3141	1589567	binsC
bin.11	69.68	11.13	0.296	Bacteria	4383	1747019	binsBC
bin.13	61.07	0.606	0.444	Bacteroidales	8236	1791288	binsAC
bin.8	54.38	7.894	0.267	Bacteria	3541	1214050	binsBC
bin.2	52.80	1.754	0.431	Bacteria	4709	2551706	binsA
bin.15	52.62	8.550	0.572	Clostridiales	1977	2116551	binsC
```

![image](https://user-images.githubusercontent.com/11639261/142347793-08b91a49-607c-4a2a-850c-3dc3ad066177.png)
- Note: *This plot compares metaWRAP performance vs. individual binning tools performance regarding both completeness and contamination*


- SRA example data (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR16888420)
  - Note: Here only one bin/MAG was recovered
```
 bin     completeness    contamination   GC      lineage N50     size    binner
 bin.1   52.44   11.26   0.532   Bacteria        1303    1111229 binsAB
```

#### 6- Run metaWRAP-REASSEMBLE_BINS
**WILL NOT BE RUN DURING DEMONSTRATION DUE TO HIGH COMPUTATIONAL COST FOR VM**

```
#metawrap reassemble_bins -o BIN_REASSEMBLY_srr -1 $CLEAN_READS/ALL_READS_1.fastq -2 $CLEAN_READS/ALL_READS_2.fastq -t 8 \
-m 300 -c 50 -x 15 -b BIN_REFINEMENT_srr/metawrap_50_15_bins

```

- Binning reassemble (metawrap example data)
```
bin	completeness	contamination	GC	lineage	N50	size
bin.10.orig	84.77	2.516	0.445	Clostridiales	2839	1452931
bin.11.permissive	70.06	10.28	0.296	Bacteria	5049	1753228
bin.12.strict	82.06	1.034	0.464	Bacteria	18867	4239208
bin.13.strict	61.65	0.652	0.444	Bacteroidales	9734	1803127
bin.15.permissive	55.27	8.885	0.572	Clostridiales	2144	2158760
bin.1.strict	71.08	1.725	0.422	Selenomonadales	3881	1658693
bin.2.orig	52.80	1.754	0.431	Bacteria	4709	2551706
bin.4.strict	70.25	0.854	0.509	Clostridiales	3450	1600391
bin.5.strict	99.32	1.006	0.406	Clostridiales	72398	2094670
bin.6.strict	85.74	1.027	0.374	Clostridiales	7014	2072694
bin.7.permissive	75.33	0.396	0.287	Clostridiales	23157	1463229
bin.8.strict	56.14	7.894	0.267	Bacteria	3859	1213933
bin.9.orig	98.4	1.6	0.311	Euryarchaeota	12258	1688675
```
![reassembled_bins](https://user-images.githubusercontent.com/11639261/142348161-a7cf60f1-09c0-4366-8440-39d2eb8ab752.png)

- Binning reassemble: SRA example data (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR16888420)
```
bin	completeness	contamination	GC	lineage	N50	size
bin.1.permissive	53.46	10.56	0.533	Bacteria	1332	1130691
```
![reassembled_bins](https://user-images.githubusercontent.com/11639261/142383636-60ee6a9b-26fd-49de-a47a-3d77dd98d5bb.png)


- **Download final bins**
``` 
#download final bins
cd $ANALYSIS_FOLDER
mkdir bins_consolidated
cd bins_consolidated

wget https://github.com/jmartinsjrbr/introduction_to_microbiomics-practice/blob/main/bin.1.permissive.fa.gz
gunzip *.gz

wget https://github.com/jmartinsjrbr/introduction_to_microbiomics-practice/blob/main/binning_demo.tar.gz
tar -vzxf *.gz
rm *.gz
cd ..
```
#### 7- Run GTDBT to assign taxonomy to bins/MAGs
```
#gtdbtk classify_wf --genome_dir bins_consolidated --out_dir gtdbtk_out --prefix omics_course.gtdbtk --cpus 10 --pplacer_cpus 1 --extension fa
```
```
user_genome     classification  fastani_reference       fastani_reference_radius        fastani_taxonomy        fastani_ani     fastani_af      closest_placement_reference     closest_placement_radius        closest_placement_taxonomy          closest_placement_ani   closest_placement_af    pplacer_taxonomy        classification_method   note    other_related_references(genome_id,species_name,radius,ANI,AF)  msa_percent     translation_table       red_value       warn    ings
bin_demo.9.orig d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter_A;s__Methanobrevibacter_A smithii GCF_000016525.1 95.0    d__Archaea;p__Methanobacteriota;c__Methano    bacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter_A;s__Methanobrevibacter_A smithii 98.04   0.92    GCF_000016525.1 95.0    d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanob    acteriaceae;g__Methanobrevibacter_A;s__Methanobrevibacter_A smithii 98.04   0.92    d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter_A;s__     taxonomic classificati    on defined by topology and ANI    topological placement and ANI have congruent species assignments        GCF_002252585.1, s__Methanobrevibacter_A smithii_A, 95.0, 92.96, 0.94; GCF_003111605.1, s__Methanobrevibacter_A woesei, 95.0, 80    .22, 0.49; GCF_900769095.1, s__Methanobrevibacter_A sp900769095, 95.0, 79.85, 0.47; GCA_902785855.1, s__Methanobrevibacter_A sp902785855, 95.0, 79.75, 0.49; GCF_003814835.1, s__Methanobrevibacter_A gottschalkii, 95.0, 79.68, 0.48; GCF    _003111625.1, s__Methanobrevibacter_A thaueri, 95.0, 79.66, 0.27; GCF_001639275.1, s__Methanobrevibacter_A oralis, 95.0, 79.63, 0.51; GCA_900319535.1, s__Methanobrevibacter_A sp900319535, 95.0, 79.58, 0.39; GCA_902774605.1, s__Methano    brevibacter_A sp902774605, 95.0, 79.54, 0.23; GCA_900314695.1, s__Methanobrevibacter_A sp900314695, 95.0, 79.42, 0.41; GCF_001548675.1, s__Methanobrevibacter_A sp001548675, 95.0, 79.4, 0.39; GCF_001477655.1, s__Methanobrevibacter_A mi    llerae_A, 95.0, 79.37, 0.46; GCA_902764455.1, s__Methanobrevibacter_A sp902764455, 95.0, 79.19, 0.39; GCA_902763935.1, s__Methanobrevibacter_A sp902763935, 95.0, 79.17, 0.34; GCF_900766745.1, s__Methanobrevibacter_A sp900766745, 95.0,     79.15, 0.34; GCA_902801725.1, s__Methanobrevibacter_A sp902801725, 95.0, 79.12, 0.36; GCA_900317865.1, s__Methanobrevibacter_A sp900317865, 95.0, 79.08, 0.42; GCA_900314615.1, s__Methanobrevibacter_A sp900314615, 95.0, 79.03, 0.39; G    CF_900103415.1, s__Methanobrevibacter_A millerae, 95.0, 79.02, 0.25; GCA_902788255.1, s__Methanobrevibacter_A sp902788255, 95.0, 78.98, 0.3; GCA_900316895.1, s__Methanobrevibacter_A sp900316895, 95.0, 78.91, 0.28; GCA_902770715.1, s__    Methanobrevibacter_A sp902770715, 95.0, 78.9, 0.35; GCA_902777885.1, s__Methanobrevibacter_A sp902777885, 95.0, 78.9, 0.31; GCA_902789505.1, s__Methanobrevibacter_A sp902789505, 95.0, 78.85, 0.34; GCA_900318035.1, s__Methanobrevibacte    r_A sp900318035, 95.0, 78.84, 0.33; GCA_902768765.1, s__Methanobrevibacter_A sp902768765, 95.0, 78.82, 0.33; GCA_902771435.1, s__Methanobrevibacter_A sp902771435, 95.0, 78.79, 0.33; GCA_900319985.1, s__Methanobrevibacter_A sp900319985    , 95.0, 78.78, 0.37; GCA_900320515.1, s__Methanobrevibacter_A sp900320515, 95.0, 78.66, 0.23; GCA_002496065.1, s__Methanobrevibacter_A sp002496065, 95.0, 78.63, 0.37; GCA_900313645.1, s__Methanobrevibacter_A sp900313645, 95.0, 78.61,     0.35; GCA_900320955.1, s__Methanobrevibacter_A sp900320955, 95.0, 78.58, 0.34; GCA_902783085.1, s__Methanobrevibacter_A sp902783085, 95.0, 78.32, 0.43; GCA_902794475.1, s__Methanobrevibacter_A sp902794475, 95.0, 77.81, 0.23; GCA_90276    4015.1, s__Methanobrevibacter_A sp902764015, 95.0, 77.73, 0.14; GCA_902796575.1, s__Methanobrevibacter_A sp902796575, 95.0, 77.63, 0.3; GCA_902772665.1, s__Methanobrevibacter_A sp902772665, 95.0, 77.26, 0.24; GCA_902782595.1, s__Metha    nobrevibacter_A sp902782595, 95.0, 76.75, 0.22    96.12   11      N/A     N/A
bin.1.permissive	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Paraclostridium;s__	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Paraclostridium;s__	ANI	N/A	GCF_006802875.1, s__Paraclostridium bifermentans, 96.2286, 99.36, 0.12; GCF_001006285.1, s__Paraclostridium benzoelyticum, 95.4733, 98.91, 0.12; GCF_000452245.2, s__Paraclostridium bifermentans_A, 96.2286, 98.64, 0.12; GCF_902362945.1, s__Paraclostridium tenue, 95.0, 89.89, 0.12; GCA_000953675.1, s__Paraclostridium sordellii, 95.0, 88.95, 0.12	30.89	11	0.9944793184213166	Genome has more than 12.5% of markers with multiple hits
bin_demo.1.strict	d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Acidaminococcales;f__Acidaminococcaceae;g__Phascolarctobacterium;s__Phascolarctobacterium faecium	GCF_003269275.1	95.0	d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Acidaminococcales;f__Acidaminococcaceae;g__Phascolarctobacterium;s__Phascolarctobacterium faecium	98.85	0.86	GCF_003269275.1	95.0	d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Acidaminococcales;f__Acidaminococcaceae;g__Phascolarctobacterium;s__Phascolarctobacterium faecium	98.85	0.86	d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Acidaminococcales;f__Acidaminococcaceae;g__Phascolarctobacterium;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_900551745.1, s__Phascolarctobacterium sp900551745, 95.0, 94.2, 0.67; GCA_900545535.1, s__Phascolarctobacterium sp900545535, 95.0, 85.36, 0.75; GCA_900544795.1, s__Phascolarctobacterium sp900544795, 95.0, 82.28, 0.63; GCA_003150755.1, s__Phascolarctobacterium sp003150755, 95.0, 78.03, 0.15; GCA_000436095.1, s__Phascolarctobacterium sp000436095, 95.0, 77.71, 0.15; GCA_900554395.1, s__Phascolarctobacterium sp900554395, 95.0, 77.41, 0.1	64.26	11	N/A	N/A
bin_demo.10.orig	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-180;s__CAG-180 sp000432435	GCA_000432435.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-180;s__CAG-180 sp000432435	98.38	0.97	GCA_000432435.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-180;s__CAG-180 sp000432435	98.38	0.97	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-180;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_002437245.1, s__CAG-180 sp002437245, 95.0, 89.39, 0.92; GCA_004556705.1, s__CAG-180 sp004556705, 95.0, 79.34, 0.2; GCA_900540295.1, s__CAG-180 sp900540295, 95.0, 77.49, 0.08; GCA_900315985.1, s__CAG-180 sp900315985, 95.0, 76.89, 0.12; GCA_004555265.1, s__CAG-180 sp004555265, 95.0, 76.74, 0.09; GCA_900545625.1, s__CAG-180 sp900545625, 95.0, 76.74, 0.08; GCA_002314305.1, s__CAG-180 sp002314305, 95.0, 76.55, 0.06; GCF_900766945.1, s__CAG-180 sp900766945, 95.0, 76.55, 0.08; GCA_002361875.1, s__CAG-180 sp002361875, 95.0, 76.07, 0.07; GCA_002490525.1, s__CAG-180 sp002490525, 95.0, 75.81, 0.08; GCA_002490425.1, s__CAG-180 sp002490425, 95.0, 75.57, 0.11	70.36	11	N/A	N/A
bin_demo.11.permissive	d__Bacteria;p__Firmicutes;c__Bacilli;o__Erysipelotrichales;f__Erysipelatoclostridiaceae;g__Faecalibacillus;s__Faecalibacillus intestinalis	GCF_003024685.1	95.0	d__Bacteria;p__Firmicutes;c__Bacilli;o__Erysipelotrichales;f__Erysipelatoclostridiaceae;g__Faecalibacillus;s__Faecalibacillus intestinalis	95.74	0.87	GCF_003024685.1	95.0	d__Bacteria;p__Firmicutes;c__Bacilli;o__Erysipelotrichales;f__Erysipelatoclostridiaceae;g__Faecalibacillus;s__Faecalibacillus intestinalis	95.74	0.87	d__Bacteria;p__Firmicutes;c__Bacilli;o__Erysipelotrichales;f__Erysipelatoclostridiaceae;g__Faecalibacillus;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_900544435.1, s__Faecalibacillus sp900544435, 95.0, 93.03, 0.72; GCF_003024675.1, s__Faecalibacillus faecis, 95.0, 90.52, 0.8; GCF_003480255.1, s__Faecalibacillus sp003480255, 95.0, 80.34, 0.52	50.29	11	N/A	Genome has more than 10.0% of markers with multiple hits
bin_demo.12.strict	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis	GCF_000154205.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis	98.19	0.82	GCF_000154205.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis	98.19	0.82	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCF_000614125.1, s__Bacteroides rodentium, 95.0, 94.39, 0.69; GCF_004793475.1, s__Bacteroides sp002491635, 95.0, 91.79, 0.63; GCA_902388495.1, s__Bacteroides sp902388495, 95.0, 81.55, 0.46; GCF_000195635.1, s__Bacteroides fluxus, 95.0, 81.47, 0.44; GCF_900108345.1, s__Bacteroides ndongoniae, 95.0, 80.76, 0.28; GCF_000155815.1, s__Bacteroides eggerthii, 95.0, 80.75, 0.37; GCF_000154525.1, s__Bacteroides stercoris, 95.0, 80.67, 0.39; GCF_900241005.1, s__Bacteroides cutis, 95.0, 80.43, 0.33; GCF_003438615.1, s__Bacteroides sp003545565, 95.0, 80.28, 0.37; GCF_900129655.1, s__Bacteroides clarus, 95.0, 80.2, 0.37; GCF_000374365.1, s__Bacteroides gallinarum, 95.0, 80.09, 0.37; GCF_004342845.1, s__Bacteroides heparinolyticus, 95.0, 79.75, 0.36; GCF_000513195.1, s__Bacteroides timonensis, 95.0, 79.68, 0.35; GCF_000186225.1, s__Bacteroides helcogenes, 95.0, 79.68, 0.37; GCA_900557355.1, s__Bacteroides sp900557355, 95.0, 79.5, 0.18; GCF_000614165.1, s__Bacteroides stercorirosoris, 95.0, 79.46, 0.36; GCF_002998435.1, s__Bacteroides zoogleoformans, 95.0, 79.35, 0.34; GCF_000315485.1, s__Bacteroides oleiciplenus, 95.0, 79.22, 0.36; GCA_900552405.1, s__Bacteroides sp900552405, 95.0, 79.08, 0.34; GCF_902364365.1, s__Bacteroides sp900556215, 95.2143, 79.07, 0.35; GCA_900555635.1, s__Bacteroides sp900555635, 95.0, 79.07, 0.39; GCF_003464595.1, s__Bacteroides intestinalis_A, 95.2143, 79.02, 0.35; GCA_900755095.1, s__Bacteroides sp900755095, 95.0, 79.02, 0.11; GCA_900765785.1, s__Bacteroides sp900765785, 95.0, 78.97, 0.17; GCF_000172175.1, s__Bacteroides intestinalis, 95.0, 78.93, 0.34; GCF_000158035.1, s__Bacteroides cellulosilyticus, 95.0, 78.77, 0.35; GCF_900155865.1, s__Bacteroides bouchesdurhonensis, 95.0, 78.67, 0.14; GCF_014334015.1, s__Bacteroides sp003463205, 95.0, 78.57, 0.19; GCF_000156195.1, s__Bacteroides finegoldii, 95.8579, 78.52, 0.19; GCF_009193295.2, s__Bacteroides sp007097645, 95.0, 78.48, 0.18; GCF_902362375.1, s__Bacteroides sp902362375, 95.0, 78.47, 0.2; GCF_000613385.1, s__Bacteroides acidifaciens, 95.0, 78.46, 0.19; GCF_001688725.2, s__Bacteroides caecimuris, 95.0, 78.42, 0.17; GCF_000011065.1, s__Bacteroides thetaiotaomicron, 95.0, 78.38, 0.2; GCA_900760255.1, s__Bacteroides sp900760255, 95.0, 78.38, 0.17; GCA_000210075.1, s__Bacteroides xylanisolvens, 95.0, 78.35, 0.18; GCF_900128475.1, s__Bacteroides massiliensis_A, 95.0, 78.27, 0.23; GCF_900130135.1, s__Bacteroides togonis, 95.0, 78.25, 0.21; GCF_009193325.2, s__Bacteroides sp900066265, 95.8579, 78.2, 0.19; GCA_900762525.1, s__Bacteroides sp900762525, 95.0, 78.19, 0.13; GCF_003865075.1, s__Bacteroides sp003865075, 95.0, 78.14, 0.14; GCF_001314995.1, s__Bacteroides ovatus, 95.0, 78.1, 0.19; GCF_002222615.2, s__Bacteroides caccae, 95.0, 78.09, 0.18; GCA_002293435.1, s__Bacteroides sp002293435, 95.0, 78.06, 0.23; GCA_900553815.1, s__Bacteroides sp900553815, 95.0, 78.06, 0.2; GCF_000025985.1, s__Bacteroides fragilis, 95.0, 78.01, 0.18; GCF_900130125.1, s__Bacteroides congonensis, 95.0, 78.0, 0.19; GCF_903181435.1, s__Bacteroides acidifaciens_A, 95.0, 77.99, 0.2; GCF_000613465.1, s__Bacteroides nordii, 95.0, 77.91, 0.17; GCF_002849695.1, s__Bacteroides fragilis_A, 95.0, 77.84, 0.18; GCF_012113595.1, s__Bacteroides sp012113595, 95.0, 77.82, 0.2; GCA_900766005.1, s__Bacteroides sp900766005, 95.0, 77.81, 0.18; GCA_900761785.1, s__Bacteroides sp900761785, 95.0, 77.72, 0.19; GCA_014385165.1, s__Bacteroides sp014385165, 95.0, 77.68, 0.15; GCF_000614145.1, s__Bacteroides faecichinchillae, 95.0, 77.67, 0.13; GCF_900106755.1, s__Bacteroides faecis, 95.0, 77.64, 0.19; GCF_000381365.1, s__Bacteroides salyersiae, 95.0, 77.6, 0.19; GCA_900547205.1, s__Bacteroides sp900547205, 95.0, 77.58, 0.17; GCA_900765805.1, s__Bacteroides sp900765805, 95.0, 77.53, 0.11; GCF_000428105.1, s__Bacteroides pyogenes, 95.0, 77.37, 0.13; GCA_900556625.1, s__Bacteroides sp900556625, 95.0, 77.34, 0.14; GCF_014196225.1, s__Bacteroides pyogenes_A, 95.0, 77.19, 0.15; GCF_000517545.1, s__Bacteroides reticulotermitis, 95.0, 77.15, 0.11; GCF_000499785.1, s__Bacteroides neonati, 95.0, 77.08, 0.13; GCF_002160055.1, s__Bacteroides sp002160055, 95.0, 77.01, 0.15; GCA_900766195.1, s__Bacteroides sp900766195, 95.0, 76.85, 0.1; GCF_007896885.1, s__Bacteroides sp007896885, 95.0, 76.83, 0.11; GCF_010500965.1, s__Bacteroides sp010500965, 95.0, 76.79, 0.11; GCA_002471185.1, s__Bacteroides sp002471185, 95.0, 76.73, 0.08; GCF_010500995.1, s__Bacteroides sp010500995, 95.0, 76.6, 0.07; GCF_010501015.1, s__Bacteroides sp010501015, 95.0, 76.56, 0.05; GCF_010500955.1, s__Bacteroides sp010500955, 95.0, 76.52, 0.07; GCA_002471195.1, s__Bacteroides sp002471195, 95.0, 76.5, 0.06; GCF_900104585.1, s__Bacteroides ihuae, 95.0, 76.46, 0.07; GCF_000428125.1, s__Bacteroides graminisolvens, 95.0, 76.46, 0.08; GCA_009929715.1, s__Bacteroides sp009929715, 95.0, 76.39, 0.09; GCA_002307035.1, s__Bacteroides sp002307035, 95.0, 76.33, 0.08	69.64	11	N/A	N/A
bin_demo.13.strict	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A	GCF_003437535.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A	97.58	0.85	GCF_003437535.1	95.0	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A	97.58	0.85	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCF_000187895.1, s__Phocaeicola plebeius, 95.0, 94.25, 0.83; GCA_900557085.1, s__Phocaeicola sp900557085, 95.0, 93.7, 0.55; GCA_900551445.1, s__Phocaeicola sp900551445, 95.0, 83.07, 0.62; GCF_902362595.1, s__Phocaeicola sp900066455, 95.0, 82.96, 0.72; GCF_900128455.1, s__Phocaeicola mediterraneensis, 95.0, 82.84, 0.7; GCA_900552075.1, s__Phocaeicola sp900552075, 95.0, 82.77, 0.62; GCA_900551065.1, s__Phocaeicola sp900551065, 95.0, 82.36, 0.65; GCF_900751965.1, s__Phocaeicola sp900751965, 95.0, 80.92, 0.18; GCF_000012825.1, s__Phocaeicola vulgatus, 95.359, 80.37, 0.25; GCF_013009555.1, s__Phocaeicola dorei, 95.359, 80.35, 0.21; GCF_900760795.1, s__Phocaeicola sp900760795, 95.0, 79.89, 0.22; GCF_000614185.1, s__Phocaeicola sartorii, 95.0, 79.4, 0.24; GCF_013618865.1, s__Phocaeicola sp900541515, 95.0, 79.36, 0.34; GCF_000374585.1, s__Phocaeicola barnesiae, 95.0, 79.14, 0.37; GCF_011959205.1, s__Phocaeicola sp011959205, 95.0, 79.07, 0.2; GCA_900551645.1, s__Phocaeicola sp900551645, 95.0, 78.99, 0.26; GCF_000157915.1, s__Phocaeicola coprophilus, 95.0, 78.88, 0.33; GCA_900553715.1, s__Phocaeicola sp900553715, 95.0, 78.8, 0.25; GCF_000154845.1, s__Phocaeicola coprocola, 95.0, 78.78, 0.29; GCA_000432735.1, s__Phocaeicola sp000432735, 95.0, 78.72, 0.34; GCA_900553185.1, s__Phocaeicola sp900553185, 95.0, 78.64, 0.29; GCF_000382445.1, s__Phocaeicola massiliensis, 95.0, 78.59, 0.23; GCA_900546645.1, s__Phocaeicola sp900546645, 95.0, 78.56, 0.28; GCF_900128495.1, s__Phocaeicola ilei, 95.0, 78.51, 0.29; GCA_900542985.1, s__Phocaeicola sp900542985, 95.0, 78.48, 0.29; GCA_900544075.1, s__Phocaeicola sp900544075, 95.0, 78.33, 0.26; GCA_900556845.1, s__Phocaeicola sp900556845, 95.0, 78.24, 0.19; GCF_902362605.1, s__Phocaeicola sp900066445, 95.0, 78.05, 0.22; GCA_902388365.1, s__Phocaeicola sp902388365, 95.0, 78.01, 0.24; GCA_000434735.1, s__Phocaeicola sp000434735, 95.0, 77.98, 0.24; GCA_000436795.1, s__Phocaeicola sp000436795, 95.0, 77.92, 0.26; GCA_002493165.1, s__Phocaeicola sp002493165, 95.0, 77.79, 0.23; GCA_900554435.1, s__Phocaeicola sp900554435, 95.0, 77.7, 0.15; GCA_900540105.1, s__Phocaeicola sp900540105, 95.0, 77.56, 0.25; GCF_000190575.1, s__Phocaeicola salanitronis, 95.0, 77.55, 0.21; GCA_900546095.1, s__Phocaeicola sp900546095, 95.0, 77.52, 0.25; GCF_002161765.1, s__Phocaeicola sp002161765, 95.0, 77.39, 0.23; GCA_900544675.1, s__Phocaeicola sp900544675, 95.0, 77.31, 0.22; GCA_900546355.1, s__Phocaeicola sp900546355, 95.0, 77.28, 0.2; GCA_004558305.1, s__Phocaeicola plebeius_B, 95.0, 77.17, 0.18; GCF_002161565.1, s__Phocaeicola sp002161565, 95.0, 77.16, 0.23; GCA_900552645.1, s__Phocaeicola sp900552645, 95.0, 76.96, 0.09; GCF_000613805.1, s__Phocaeicola paurosaccharolyticus, 95.0, 76.87, 0.07; GCA_902779535.1, s__Phocaeicola sp902779535, 95.0, 76.25, 0.06; GCF_000312445.1, s__Phocaeicola abscessus, 95.0, 76.03, 0.04; GCA_002315285.1, s__Phocaeicola sp002315285, 95.0, 75.51, 0.03	47.39	11	N/A	N/A
bin_demo.15.permissive	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__ER4;s__	N/A	N/A	N/A	N/A	N/A	GCF_000765235.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__ER4;s__ER4 sp000765235	95.85	0.52	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__ER4;s__	taxonomic classification defined by topology and ANI	N/A	GCA_900552015.1, s__ER4 sp900552015, 95.0, 89.28, 0.5; GCA_900550165.1, s__ER4 sp900550165, 95.0, 87.74, 0.45; GCA_003522105.1, s__ER4 sp003522105, 95.0, 87.01, 0.46; GCA_900763705.1, s__ER4 sp900763705, 95.0, 86.75, 0.35; GCA_900546295.1, s__ER4 sp900546295, 95.0, 85.93, 0.46; GCA_900556145.1, s__ER4 sp900556145, 95.0, 84.62, 0.39; GCA_002405995.1, s__ER4 sp002405995, 95.0, 83.53, 0.41; GCA_900768135.1, s__ER4 sp900768135, 95.0, 82.9, 0.39; GCA_900317525.1, s__ER4 sp900317525, 95.0, 79.81, 0.36; GCA_002437735.1, s__ER4 sp002437735, 95.0, 79.79, 0.28; GCA_900552375.1, s__ER4 sp900552375, 95.0, 79.13, 0.28	33.59	11	0.9976355171122855	N/A
bin_demo.2.orig	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter faecis	GCF_001406815.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter faecis	98.35	0.82	GCA_900557055.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter sp900557055	97.15	0.66	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__	ANI	topological placement and ANI have incongruent species assignments	GCA_002474415.1, s__Agathobacter sp002474415, 95.0, 90.5, 0.59; GCA_900550845.1, s__Agathobacter sp900550845, 95.0, 86.05, 0.17; GCF_902363675.1, s__Agathobacter sp000434275, 95.0, 81.35, 0.4; GCA_900550545.1, s__Agathobacter sp900550545, 95.0, 80.31, 0.32; GCA_900552085.1, s__Agathobacter sp900552085, 95.0, 79.56, 0.28; GCF_000020605.1, s__Agathobacter rectalis, 95.0, 79.36, 0.15; GCA_900543445.1, s__Agathobacter sp900543445, 95.0, 78.47, 0.24; GCA_900549895.1, s__Agathobacter sp900549895, 95.0, 78.29, 0.24; GCA_900547695.1, s__Agathobacter sp900547695, 95.0, 78.18, 0.13; GCA_900317585.1, s__Agathobacter sp900317585, 95.0, 78.15, 0.1; GCA_900546625.1, s__Agathobacter sp900546625, 95.0, 77.55, 0.14; GCA_900316805.1, s__Agathobacter sp900316805, 95.0, 77.16, 0.19; GCA_900548765.1, s__Agathobacter sp900548765, 95.0, 76.24, 0.07; GCF_002735305.1, s__Agathobacter ruminis, 95.0, 76.13, 0.07	42.39	11	N/A	N/A
bin_demo.4.strict	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-177;s__CAG-177 sp003514385	GCA_003514385.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-177;s__CAG-177 sp003514385	95.6	0.86	N/A	N/A	N/A	N/A	N/A	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-177;s__	ANI	topological placement and ANI have incongruent species assignments	GCA_003514385.1, s__CAG-177 sp003514385, 95.0, 95.6, 0.86; GCA_002451755.1, s__CAG-177 sp002451755, 95.0, 94.09, 0.72; GCA_002438685.1, s__CAG-177 sp002438685, 95.0, 93.15, 0.84; GCA_003538135.1, s__CAG-177 sp003538135, 95.0, 87.88, 0.85; GCA_000431775.1, s__CAG-177 sp000431775, 95.0, 82.89, 0.21; GCF_900770255.1, s__CAG-177 sp900770255, 95.0, 79.82, 0.56; GCF_900771185.1, s__CAG-177 sp900771185, 95.0, 78.14, 0.17; GCA_002309055.1, s__CAG-177 sp002309055, 95.0, 77.11, 0.02; GCA_002315935.1, s__CAG-177 sp002315935, 95.0, 75.89, 0.05; GCA_900319455.1, s__CAG-177 sp900319455, 95.0, 75.57, 0.02	48.76	11	N/A	N/A
bin_demo.5.strict	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B	GCF_002834235.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B	97.93	0.85	GCF_002834235.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B	97.93	0.85	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCF_002834225.1, s__Ruminococcus_E bromii, 95.0, 94.9, 0.82; GCA_902771755.1, s__Ruminococcus_E sp902771755, 95.0, 92.48, 0.8; GCA_902766775.1, s__Ruminococcus_E sp902766775, 95.0, 90.54, 0.81; GCF_003438075.1, s__Ruminococcus_E sp003438075, 95.0, 80.16, 0.41; GCF_014288065.1, s__Ruminococcus_E sp003526955, 95.0, 79.92, 0.37; GCA_002491825.1, s__Ruminococcus_E sp002491825, 95.0, 79.88, 0.51; GCA_003521625.1, s__Ruminococcus_E sp003521625, 95.0, 79.55, 0.45; GCA_902766885.1, s__Ruminococcus_E sp902766885, 95.0, 79.46, 0.44; GCA_002493005.1, s__Ruminococcus_E sp002493005, 95.0, 79.45, 0.43; GCF_900101355.1, s__Ruminococcus_E bromii_A, 95.0, 79.34, 0.32; GCF_900755995.1, s__Ruminococcus_E sp900755995, 95.0, 79.1, 0.29; GCA_002493635.1, s__Ruminococcus_E sp002493635, 95.0, 79.03, 0.35; GCA_900543095.1, s__Ruminococcus_E sp900543095, 95.0, 78.48, 0.36; GCA_900315085.1, s__Ruminococcus_E sp900315085, 95.0, 78.21, 0.22; GCF_900100595.1, s__Ruminococcus_E sp900100595, 95.0, 77.74, 0.11; GCA_900315605.1, s__Ruminococcus_E sp900315605, 95.0, 77.72, 0.18; GCA_900318495.1, s__Ruminococcus_E sp900318495, 95.0, 77.64, 0.23; GCA_900316555.1, s__Ruminococcus_E sp900316555, 95.0, 77.61, 0.15; GCA_900315785.1, s__Ruminococcus_E sp900315785, 95.0, 77.58, 0.13; GCA_900316435.1, s__Ruminococcus_E sp900316435, 95.0, 77.58, 0.15; GCA_902798605.1, s__Ruminococcus_E sp902798605, 95.0, 77.57, 0.14; GCA_902767845.1, s__Ruminococcus_E sp902767845, 95.0, 77.56, 0.11; GCA_004560275.1, s__Ruminococcus_E sp004560275, 95.0, 77.48, 0.19; GCA_900315195.1, s__Ruminococcus_E sp900315195, 95.0, 77.45, 0.19; GCA_902794395.1, s__Ruminococcus_E sp902794395, 95.0, 77.44, 0.09; GCA_900319655.1, s__Ruminococcus_E sp900319655, 95.0, 77.39, 0.25; GCF_005601135.1, s__Ruminococcus_E sp900314705, 95.0, 77.36, 0.12; GCA_902774325.1, s__Ruminococcus_E sp902774325, 95.0, 77.33, 0.08; GCA_900313895.1, s__Ruminococcus_E sp900313895, 95.0, 77.29, 0.02; GCA_902776375.1, s__Ruminococcus_E sp902776375, 95.0, 77.27, 0.1; GCA_902792105.1, s__Ruminococcus_E sp902792105, 95.0, 77.24, 0.18; GCA_003520555.1, s__Ruminococcus_E sp003520555, 95.0, 77.23, 0.2; GCA_902798365.1, s__Ruminococcus_E sp902798365, 95.0, 77.22, 0.13; GCA_902772755.1, s__Ruminococcus_E sp902772755, 95.0, 77.19, 0.09; GCA_902785365.1, s__Ruminococcus_E sp902785365, 95.0, 77.18, 0.11; GCA_902799325.1, s__Ruminococcus_E sp902799325, 95.0, 77.13, 0.11; GCA_002350765.1, s__Ruminococcus_E sp002350765, 95.0, 77.12, 0.09; GCA_902790475.1, s__Ruminococcus_E sp902790475, 95.0, 77.1, 0.12; GCA_902784855.1, s__Ruminococcus_E sp902784855, 95.0, 77.07, 0.08; GCA_902776975.1, s__Ruminococcus_E sp902776975, 95.0, 77.04, 0.03; GCA_900316385.1, s__Ruminococcus_E sp900316385, 95.0, 77.03, 0.11; GCA_900314795.1, s__Ruminococcus_E sp900314795, 95.0, 76.99, 0.04; GCA_900318905.1, s__Ruminococcus_E sp900318905, 95.0, 76.98, 0.17; GCA_902777275.1, s__Ruminococcus_E sp902777275, 95.0, 76.97, 0.13; GCA_902793125.1, s__Ruminococcus_E sp902793125, 95.0, 76.94, 0.1; GCA_902787145.1, s__Ruminococcus_E sp902787145, 95.0, 76.93, 0.08; GCA_900316815.1, s__Ruminococcus_E sp900316815, 95.0, 76.76, 0.06; GCA_902765625.1, s__Ruminococcus_E sp902765625, 95.0, 76.62, 0.1; GCA_902797225.1, s__Ruminococcus_E sp902797225, 95.0, 76.6, 0.07; GCA_902797655.1, s__Ruminococcus_E sp902797655, 95.0, 76.59, 0.11; GCA_900319615.1, s__Ruminococcus_E sp900319615, 95.0, 76.52, 0.03; GCA_902761375.1, s__Ruminococcus_E sp902761375, 95.0, 76.51, 0.09; GCA_900317595.1, s__Ruminococcus_E sp900317595, 95.0, 76.49, 0.09; GCA_900317875.1, s__Ruminococcus_E sp900317875, 95.0, 76.44, 0.13; GCA_002353935.1, s__Ruminococcus_E sp002353935, 95.0, 76.38, 0.06; GCA_902797525.1, s__Ruminococcus_E sp902797525, 95.0, 76.35, 0.13; GCA_902785355.1, s__Ruminococcus_E sp902785355, 95.0, 76.34, 0.04; GCA_902773635.1, s__Ruminococcus_E sp902773635, 95.0, 76.33, 0.03; GCA_900317315.1, s__Ruminococcus_E sp900317315, 95.0, 76.3, 0.08; GCA_900320415.1, s__Ruminococcus_E sp900320415, 95.0, 76.18, 0.08; GCA_902773155.1, s__Ruminococcus_E sp902773155, 95.0, 76.15, 0.04; GCA_902767325.1, s__Ruminococcus_E sp902767325, 95.0, 76.1, 0.14; GCA_902764715.1, s__Ruminococcus_E sp902764715, 95.0, 75.97, 0.03; GCA_902773035.1, s__Ruminococcus_E sp902773035, 95.0, 75.9, 0.08; GCA_902790985.1, s__Ruminococcus_E sp902790985, 95.0, 75.85, 0.03; GCA_902763855.1, s__Ruminococcus_E sp902763855, 95.0, 75.79, 0.03; GCA_900315815.1, s__Ruminococcus_E sp900315815, 95.0, 75.76, 0.03; GCA_902776185.1, s__Ruminococcus_E sp902776185, 95.0, 75.66, 0.03; GCA_002394725.1, s__Ruminococcus_E sp002394725, 95.0, 75.65, 0.03; GCA_902798755.1, s__Ruminococcus_E sp902798755, 95.0, 75.51, 0.03; GCA_003499325.1, s__Ruminococcus_E sp003499325, 95.0, 75.37, 0.03; GCA_902779225.1, s__Ruminococcus_E sp902779225, 95.0, 75.25, 0.03; GCA_900320995.1, s__Ruminococcus_E sp900320995, 95.0, 75.24, 0.02; GCA_902796225.1, s__Ruminococcus_E sp902796225, 95.0, 74.97, 0.01; GCA_902797825.1, s__Ruminococcus_E sp900313865, 95.0, 74.89, 0.01	96.64	11	N/A	N/A
bin_demo.6.strict	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Butyrivibrio_A;s__Butyrivibrio_A crossotus	GCF_000156015.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Butyrivibrio_A;s__Butyrivibrio_A crossotus	98.65	0.88	GCF_000156015.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Butyrivibrio_A;s__Butyrivibrio_A crossotus	98.65	0.88	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Butyrivibrio_A;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_900543865.1, s__Butyrivibrio_A sp900543865, 95.0, 98.04, 0.7; GCA_000431815.1, s__Butyrivibrio_A sp000431815, 95.0, 77.71, 0.19; GCF_900768755.1, s__Butyrivibrio_A sp900768755, 95.0, 76.61, 0.08; GCF_900771195.1, s__Butyrivibrio_A sp900771195, 95.0, 76.12, 0.07	65.61	11	N/A	N/A
bin_demo.7.permissive	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-508;g__CAG-273;s__CAG-273 sp000435755	GCA_000435755.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-508;g__CAG-273;s__CAG-273 sp000435755	99.92	0.89	GCA_000435755.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-508;g__CAG-273;s__CAG-273 sp000435755	99.92	0.89	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-508;g__CAG-273;s__	taxonomic classification defined by topology and ANI	topological placement and ANI have congruent species assignments	GCA_000437855.1, s__CAG-273 sp000437855, 95.0, 86.98, 0.65; GCF_900752095.1, s__CAG-273 sp900752095, 95.0, 79.46, 0.35; GCF_900752335.1, s__CAG-273 sp900752335, 95.0, 79.12, 0.29; GCA_902406115.1, s__CAG-273 sp902406115, 95.0, 78.36, 0.35; GCA_003507395.1, s__CAG-273 sp003507395, 95.0, 78.06, 0.35; GCA_003534295.1, s__CAG-273 sp003534295, 95.0, 77.99, 0.29; GCA_000438355.1, s__CAG-273 sp000438355, 95.0, 77.92, 0.33	87.51	11	N/A	N/A
bin_demo.8.strict	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-465;g__CAG-465;s__	N/A	N/A	N/A	N/A	N/A	GCA_000433755.1	95.0	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-465;g__CAG-465;s__CAG-465 sp000433755	97.15	0.63	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__TANB77;f__CAG-465;g__CAG-465;s__	taxonomic classification defined by topology and ANI	N/A	GCF_900772475.1, s__CAG-465 sp900772475, 95.0, 93.48, 0.62; GCA_900554875.1, s__CAG-465 sp900554875, 95.0, 88.46, 0.53; GCA_000433335.1, s__CAG-465 sp000433335, 95.0, 77.1, 0.18; GCA_000433135.1, s__CAG-465 sp000433135, 95.0, 76.9, 0.21	46.42	11	0.997839921360666	N/A

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

- GFF file resulted from prokka
```
NODE_314_length_15207_cov_0.831639	prokka	gene	168	1130	.	-	.	ID=JIEOEPFJ_00001_gene;Name=pfkA;gene=pfkA;locus_tag=JIEOEPFJ_00001
NODE_314_length_15207_cov_0.831639	prokka	mRNA	168	1130	.	-	.	ID=JIEOEPFJ_00001_mRNA;Name=pfkA;gene=pfkA;locus_tag=JIEOEPFJ_00001
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	168	1130	.	-	0	ID=JIEOEPFJ_00001;Parent=JIEOEPFJ_00001_gene,JIEOEPFJ_00001_mRNA;eC_number=2.7.1.11;Name=pfkA;gene=pfkA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P00512;locus_tag=JIEOEPFJ_00001;product=ATP-dependent 6-phosphofructokinase
NODE_314_length_15207_cov_0.831639	prokka	gene	1212	4718	.	-	.	ID=JIEOEPFJ_00002_gene;Name=dnaE;gene=dnaE;locus_tag=JIEOEPFJ_00002
NODE_314_length_15207_cov_0.831639	prokka	mRNA	1212	4718	.	-	.	ID=JIEOEPFJ_00002_mRNA;Name=dnaE;gene=dnaE;locus_tag=JIEOEPFJ_00002
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	1212	4718	.	-	0	ID=JIEOEPFJ_00002;Parent=JIEOEPFJ_00002_gene,JIEOEPFJ_00002_mRNA;eC_number=2.7.7.7;Name=dnaE;db_xref=COG:COG0587;gene=dnaE;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:Q9XDH5;locus_tag=JIEOEPFJ_00002;product=DNA polymerase III subunit alpha
NODE_314_length_15207_cov_0.831639	prokka	gene	4735	6033	.	-	.	ID=JIEOEPFJ_00003_gene;Name=lysA;gene=lysA;locus_tag=JIEOEPFJ_00003
NODE_314_length_15207_cov_0.831639	prokka	mRNA	4735	6033	.	-	.	ID=JIEOEPFJ_00003_mRNA;Name=lysA;gene=lysA;locus_tag=JIEOEPFJ_00003
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	4735	6033	.	-	0	ID=JIEOEPFJ_00003;Parent=JIEOEPFJ_00003_gene,JIEOEPFJ_00003_mRNA;eC_number=4.1.1.20;Name=lysA;db_xref=COG:COG0019;gene=lysA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P9WIU7;locus_tag=JIEOEPFJ_00003;product=Diaminopimelate decarboxylase
NODE_314_length_15207_cov_0.831639	prokka	gene	6048	6626	.	-	.	ID=JIEOEPFJ_00004_gene;Name=engB;gene=engB;locus_tag=JIEOEPFJ_00004
NODE_314_length_15207_cov_0.831639	prokka	mRNA	6048	6626	.	-	.	ID=JIEOEPFJ_00004_mRNA;Name=engB;gene=engB;locus_tag=JIEOEPFJ_00004
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	6048	6626	.	-	0	ID=JIEOEPFJ_00004;Parent=JIEOEPFJ_00004_gene,JIEOEPFJ_00004_mRNA;Name=engB;db_xref=COG:COG0218;gene=engB;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P38424;locus_tag=JIEOEPFJ_00004;product=putative GTP-binding protein EngB
NODE_314_length_15207_cov_0.831639	prokka	gene	6637	9054	.	-	.	ID=JIEOEPFJ_00005_gene;Name=lon1;gene=lon1;locus_tag=JIEOEPFJ_00005
NODE_314_length_15207_cov_0.831639	prokka	mRNA	6637	9054	.	-	.	ID=JIEOEPFJ_00005_mRNA;Name=lon1;gene=lon1;locus_tag=JIEOEPFJ_00005
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	6637	9054	.	-	0	ID=JIEOEPFJ_00005;Parent=JIEOEPFJ_00005_gene,JIEOEPFJ_00005_mRNA;eC_number=3.4.21.53;Name=lon1;db_xref=COG:COG0466;gene=lon1;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P37945;locus_tag=JIEOEPFJ_00005;product=Lon protease 1
NODE_314_length_15207_cov_0.831639	prokka	gene	9302	10384	.	+	.	ID=JIEOEPFJ_00006_gene;Name=serC;gene=serC;locus_tag=JIEOEPFJ_00006
NODE_314_length_15207_cov_0.831639	prokka	mRNA	9302	10384	.	+	.	ID=JIEOEPFJ_00006_mRNA;Name=serC;gene=serC;locus_tag=JIEOEPFJ_00006
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	9302	10384	.	+	0	ID=JIEOEPFJ_00006;Parent=JIEOEPFJ_00006_gene,JIEOEPFJ_00006_mRNA;eC_number=2.6.1.52;Name=serC;db_xref=COG:COG1932;gene=serC;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:Q9HZ66;locus_tag=JIEOEPFJ_00006;product=Phosphoserine aminotransferase
NODE_314_length_15207_cov_0.831639	prokka	gene	10402	11565	.	+	.	ID=JIEOEPFJ_00007_gene;Name=serA;gene=serA;locus_tag=JIEOEPFJ_00007
NODE_314_length_15207_cov_0.831639	prokka	mRNA	10402	11565	.	+	.	ID=JIEOEPFJ_00007_mRNA;Name=serA;gene=serA;locus_tag=JIEOEPFJ_00007
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	10402	11565	.	+	0	ID=JIEOEPFJ_00007;Parent=JIEOEPFJ_00007_gene,JIEOEPFJ_00007_mRNA;eC_number=1.1.1.95;Name=serA;db_xref=COG:COG0111;gene=serA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P0A9T0;locus_tag=JIEOEPFJ_00007;product=D-3-phosphoglycerate dehydrogenase
NODE_314_length_15207_cov_0.831639	prokka	gene	11643	12518	.	-	.	ID=JIEOEPFJ_00008_gene;Name=polC;gene=polC;locus_tag=JIEOEPFJ_00008
NODE_314_length_15207_cov_0.831639	prokka	mRNA	11643	12518	.	-	.	ID=JIEOEPFJ_00008_mRNA;Name=polC;gene=polC;locus_tag=JIEOEPFJ_00008
NODE_314_length_15207_cov_0.831639	Prodigal:002006	CDS	11643	12518	.	-	0	ID=JIEOEPFJ_00008;Parent=JIEOEPFJ_00008_gene,JIEOEPFJ_00008_mRNA;eC_number=2.7.7.7;Name=polC;gene=polC;inference=ab initio prediction:Prodigal
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
- run_dbcan results for CAZymes prediction by running hmmer algorithm
```
GeneID  HMMER Hotpep  DIAMOND #ofTools
DEBFJKAO_00183	GT66(22-798)	-	-	1
DEBFJKAO_00449	GT4(188-336)	-	-	1
DEBFJKAO_00524	GT2_Glycos_transf_2(5-159)	-	-	1
DEBFJKAO_00527	GT81(19-243)	-	-	1
DEBFJKAO_00531	GT2_Glycos_transf_2(12-150)	-	-	1
DEBFJKAO_00695	GT2_Glycos_transf_2(6-137)	-	-	1
DEBFJKAO_00699	GT2_Glycos_transf_2(5-124)	-	-	1
DEBFJKAO_00753	AA3(2-398)	-	-	1
DEBFJKAO_00956	GT2_Glycos_transf_2(6-127)	-	-	1
DEBFJKAO_00959	GT2_Glycos_transf_2(139-270)	-	-	1
DEBFJKAO_00960	GT2_Glycos_transf_2(6-130)	-	-	1
DEBFJKAO_00961	GT2_Glycos_transf_2(5-146)	-	-	1
DEBFJKAO_01087	GH109(128-271)	-	-	1
DEBFJKAO_01112	GT2_Glycos_transf_2(4-111)	-	-	1
DEBFJKAO_01121	GT2_Glycos_transf_2(6-139)	-	-	1
DEBFJKAO_01165	GH109(2-149)	-	-	1
DEBFJKAO_01277	CE1(11-261)	-	-	1
DEBFJKAO_01363	GT2_Glycos_transf_2(7-178)+GT2_Glycos_transf_2(338-476)	-	-	1
DEBFJKAO_01437	GT2_Glycos_transf_2(4-129)	-	-	1
DEBFJKAO_01438	GT2_Glycos_transf_2(4-178)	-	-	1
DEBFJKAO_01439	GT2_Glycos_transf_2(7-108)	-	-	1
DEBFJKAO_01468	GT2_Glycos_transf_2(5-82)	-	-	1
DEBFJKAO_01473	GT2_Glycos_transf_2(19-179)	-	-	1
DEBFJKAO_01532	GT2_Glycos_transf_2(136-257)	-	-	1
DEBFJKAO_01534	GT2_Glycos_transf_2(6-129)	-	-	1
DEBFJKAO_01554	GH154(15-368)	-	-	1
DEBFJKAO_01555	GH154(27-379)	-	-	1
```
