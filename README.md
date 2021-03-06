DPAG fly TRAP 
Kevin Talbot 'Shatra'
Jakob Sca
  - See minion repo

# This repo is now deprecated and is replaced by version 2 of the analysis, found [here](https://github.com/peterkilfeather/trap).

## 15th May 2020
Ongoing routes:
  ### Deconvolution of axonal TRAP data:
  - After reading a [benchmarking comparison](https://www.biorxiv.org/content/10.1101/2020.01.10.897116v1.full) I have prioritised testing CIBERSORTx to deconvolute the axonal IP data. I have also started testing 'nnls', which I believe works in a similar method. Another tool to try is FARDEEP. I started by using RNA Seq data from Chen, 2014 (GSE52564), and found that axonal trap samples resemble 70% neurons, and the remainder astrocytes and newly-formed oligodendrocytes. I need to test what proportion NNLS and a tool like FARDEEP estimates with this Chen dataset. The Chen dataset was run through the same kallisto settings as the TRAP data. I then tested a single cell dataset from [Gokce, 2017](https://www.ncbi.nlm.nih.gov/pubmed/27425622) and it predicted 90-95% neuronal, remainder astrocyte/oligo.
  - I now want to use TRAP datasets of different neuronal types + astrocytes and oligodendrocytes to further tease apart the neuronal category. I am using the following datasets:
    - Astrocytes: [Clarke, 2018](https://www.pnas.org/content/115/8/E1896)
    - Astrocytes: [Sakers, 2017](https://www.pnas.org/content/114/19/E3830)
    - MSNs: [Kronman, 2018](https://www.biorxiv.org/content/10.1101/444315v1)
    - DAns: [Brichta, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4763340/)
    - Oligos: [Voskuhl, 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118451)
  - I plan on importing them all in as one DESeq2 object, export the counts and feed them into CIBERSORTx. **Should I import the axonal data into the same DESeq object?**
### DEG of ageing in midbrain
  - I need to test ageing datasets for differential expression and find the overlap of their genes with our 171.
    - First dataset to test: Astrocytes: [Clarke, 2018](https://www.pnas.org/content/115/8/E1896)

### Somatic variant calling in midbrain

  - I am following the off-label tutorial for calling somatic mutation differences between to samples from [GATK](https://gatkforums.broadinstitute.org/gatk/discussion/11315/off-label-workflow-to-simply-call-differences-in-two-samples).
  - So far, I have called each KW TRAP sample in tumor-only mode with Mutect2, created a panel of normals from these samples, created modified mutect2 vcf files by moving the allele frequency/fraction and have run haplotypeCaller (with sentieon) on all these samples to identify germline variants to exclude from the analysis.
  - I need to filter the output of haplotypeCaller: Can I use VQSR for this or do I need to do just hard filtering? I want to export the raw haplotypeCaller calls into R to plot the QC distributions like in [this post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).
  - I need to decide on what I am asking and how to structure the comparisons:
    - Does midbrain tissue acquire more somatic mutations with age
    - Do dopaminergic neurons have more somatic mutations, regardless of age
    - Do dopaminergic neurons aquire more somatic mutations with age
    - Do dopaminergic neurons aquire more somatic mutations with OVX, regardless of age
    - Do dopaminergic neurons acquire more somatic mutations with OVX, with age
- Splicing in midbrain

## GATK
```
remove alternate contigs:
grep -v "random" /home/peter/may_2020/mm10_protCod_canonical.bed | grep -v "JH584" > /home/peter/may_2020/mm10_protCod_canonical_noAlt.bed

Get all INFO and FORMAT fields from a VCF:
zcat 257_273_haplo_gatk_hardfilt.vcf.gz | perl -ne '/^.*INFO.*ID=([a-zA-Z0-9_]+),/ && print "-F $1 \\ "' | uniq
zcat 257_273_haplo_gatk_hardfilt.vcf.gz | perl -ne '/^.*FORMAT.*ID=([a-zA-Z0-9_]+),/ && print "-GF $1 \\ "' | uniq
```

## Leafcutter
```
for i in /zfs/analysis/pk_trap/bam_april2020/*sorted.bam ; do ./scripts/bam2junc.sh $i ${i%.bam}.junc ; done
python /home/peter/april_2020/trap/leafcutter/clustering/leafcutter_cluster.py -j /zfs/analysis/pk_trap/bam_april2020/junc_files.txt -m 50 -o /zfs/analysis/pk_trap/bam_april2020/intron_cluster -l 500000


```
 
## 22nd - 28th April 2020
Examining difference between young and old midbrain samples. Plan to look at DGE, DTU/Splicing, Mutation level, 3' UTR usage and length:
  - DGE example: [Normal aging induces A1-like astrocyte reactivity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5828643/)
  - Splicing example: [Integrative transcriptome analyses of the aging brain implicate altered splicing in Alzheimer’s disease susceptibility](https://www.nature.com/articles/s41588-018-0238-1)
  - Mutation example: [Aging and neurodegeneration are associated with increased mutations in single human neurons](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5831169/)
  
  
```
for i in *1.fastq.gz ; do fq1=$i ; fq2=${i%1.fastq.gz}2.fastq.gz ; sample_code=$(echo $i | cut -d '_' -f 3) ; output=${i%_1.fastq.gz}.bam ; zcat $i | head -n 1 | awk -v sample_code=$sample_code -v fq1=$fq1 -v fq2=$fq2 -v output=$output -F":" '{ID = $3 "." $4} {PU = $3 "." $4 "." sample_code} {SM = sample_code} {PL = "ILLUMINA"} {LB = sample_code ".LIB1"} {OFS=""} { print "picard FastqToSam FASTQ=" fq1 " FASTQ2=" fq2 " OUTPUT=" output " READ_GROUP_NAME=" ID " SAMPLE_NAME=" sample_code " LIBRARY_NAME=" LB " PLATFORM_UNIT=" PU " PLATFORM=" PL " SEQUENCING_CENTER=WTCHG" }'; done > fastq_to_sam.sh
```

```
for i in *.bam ; do echo STAR   --runThreadN 8   --genomeDir /home/peter/april_2020/star_references/gencode_vM24_first_pass_PK1/   --genomeLoad LoadAndKeep   --readFilesType SAM   PE      --readFilesIn $i      --readFilesCommand samtools   view      --outFileNamePrefix /home/peter/april_2020/star/${i%.bam}_pass_2_   --outReadsUnmapped Fastx   --outSAMtype BAM   Unsorted      --outSAMstrandField intronMotif   --sjdbOverhang 74   --twopassMode None --outSAMattrRGline $(samtools view -H $i | grep @RG | awk '{OFS=" "} {print $2,$3,$4,$5,$6}') ; done
```
  
```
sed '/^\(GL\|JH\)/d' first_pass_SJ_out_KW.tab | sed 's/^/chr/' > first_pass_SJ_out_KW.tab2
```
```
for i in $(ls *.bam | cut -d _ -f 3 | cut -d . -f 1 | sort | uniq ) ; do bam_input="" ; for j in $(ls *${i}.bam) ; do bam_input="$bam_input I=$j" ; done ; echo picard MarkDuplicates $bam_input O=/tank/${i}_dedup.bam M=/tank/${i}_dup_metrics.txt ; done > mark_dups_command_list

#then feed into GNU parallel
```

```
cat P170742_B01\&B02_SampleSheet.csv | tr -dc '[:print:]\n' | awk  'NR>11 {print $0}' | awk 'BEGIN{FS=",";OFS=""}{print $1, "\t", $1}' | tr [a-z] [A-Z] | sed '1i Sample ID\tsample_ID' > /zfs/analysis/neurochip/sample_code_to_id_section1.txt
cat P170742_B03\&B04_SampleSheet.csv | tr -dc '[:print:]\n' | awk  'NR>11 {print $0}' | awk 'BEGIN{FS=",";OFS=""}{print $4, "_", $5, "\t", $1}' | tr [a-z] [A-Z] > /zfs/analysis/neurochip/sample_code_to_id_section2.txt
cat P170742_B05\&B06_SampleSheet.csv | tr -dc '[:print:]\n' | awk  'NR>11 {print $0}' | awk 'BEGIN{FS=",";OFS=""}{print $4, "_", $5, "\t", $1}' | tr [a-z] [A-Z] > /zfs/analysis/neurochip/sample_code_to_id_section3.txt
cat sample_code_to_id_section1.txt sample_code_to_id_section2.txt sample_code_to_id_section3.txt > sample_code_to_id.txt
```

![](https://github.com/peterkilfeather/pk_trap/blob/master/TRAP.png)
 
## 16th March 2020
See 11th March for updates.

- **LCM/ST**: Over weekend, started optimising protocol for rapid DAB staining of TH, following [LCM-Seq protocol](https://www.ncbi.nlm.nih.gov/pubmed/29130192), published in [Nature comms, 2016](https://www.nature.com/articles/ncomms12139). 
  - 25 second cresyl violet exposure + 1/25 primary and secondary antibody concentrations work, 4 min incubation with antibodies and ABC solution, DAB. 30 second PBS washes, instead of 2 min.
  - Imaged on EVOS: dopamine neurons clearly labelled
  - Extracted RNA from sections and obtained 400 ng from 4 sections total.
    - RIN from fresh section: 8.8
    - RIN after staining: 5.3
  - Next step to use 1U/ul rRNasin in antibodies, ABC and DAB solution, plus use cold PBS for washes.
  - Ordered polyvinylsulfonic acid to test as cheap RNase inhibitor from [Earl, 2018 paper](https://www.ncbi.nlm.nih.gov/pubmed/28662363). 
  - Also, could investigate shorter wash conditions (dip in PBS tubes 1 and 2, 30 seconds in tube 3) and shorter incubation conditions with antibodies. 
  - Also, could test HRP-conjugated secondary.
 
## 12th March 2020
- LCM staining protocol
  - Start point with sections on slides.
    ```
    - Prepare 15 ml tubes
    - Prepare 3 x falcons with PBS + 1 x ddH2O
    - Prepare 1 x falcon with ice-cold acetone
    - Prepare 0.25% Triton X-100 in PBS (500 ul per slide)
    - Prepare primary antibody
    - Prepare secondary antibody
    - Prepare 2.5 ml of ABC solution (covered with foil and in motion for 30 min before use)
    - Prepare 2.5 ml of DAB solution (covered in foil and in motion for 2-3 min before use)
    - Dissolve solid cresyl violet acetate (e.g. ALDRICH #86,098-0) at a concentration of 1% (w/v) in 75% EtOH at room temperature with agitation/stirring for several hours to overnight
    - Adjust pH of cresyl violet solution to pH 8.0 with 3M Tris-HCl
    - Filter the staining solution before use to remove unsolubilized powder
    - (Sometimes Lot to Lot variations in the purchased cresyl violet powder can lead to weaker staining results if the dye content is below 75%)
    
    1.  Place slide in ice cold acetone for 5 min
    2.  Dry slide on lint-free tissue after every solution step
    3.  Wash slide in PBS 1 min, 3 x
    4.  Incubate ~250 ul of rabbit anti-TH 1/25 in PBS with 0.25% Triton X-100 for 4 min
    5.  Wash slide in PBS 2 min, 3 x
    6.  Incubate ~250 ul of biotinylated anti-rabbit 1/25 in PBS with 0.25% Triton X-100 for 4 min
    7.  Wash slide in PBS 2 min, 3 x
    8.  Incubate ~250 ul of ABC solution for 4 min
    9.  Wash slide in PBS 2 min, 3 x
    10. Incubate ~250 ul of DAB solution for 1-2 min
    11. Wash slide in PBS 2 min, 2 x
    12. Incubate ~250 ul of cresyl violet solution for 20-60 s
    13. Wash slide in ddH2O for 30 s
    14. Dehydrate by moving slide through 50, 75, 95 and 99.7% EtOH, 30 s each
    15. Air-dry slides for 3 min
    16. Proceed to LCM
    ```
## 11th March 2020
- Meeting with RWM on 16th at 9:30
  - Will discuss images of striatum GFP taken by GH. *16th: Discussed: RWM wants samples to be collected and frozen. Investigate strategies to try to deconvolute dopaminergic signal from the striatal TRAP signal. So looking at doing at least 8x TRAPs with all three regions, spread over 4 weeks. Currently have 3x working matrix made, require 10x matrix for 10x three-region TRAPs using 0.125x MB and 0.0625x DS/VS. Plan to make 8x matrix Tuesday morning and test on 4 REL mice. Keep the samples, as will prove useful in validation/optimisation of long-read/polyA protocols. Also freezing 5 mice brains in OCT.*
  - Priorities:
    - Start putting together powerpoint and figures of ELISA and TRAPs from this year
      - Perform second ELISA on first batch of samples
        - Running overnight *16th: Done and data added to R project for plotting*
    - Sample collection:
      - TRAP on hold until after RWM meeting
      - Need to reallocate mice to frozen (OCT) samples *Done, see [cohort design](https://docs.google.com/spreadsheets/d/1G7gFk-47ZYPr35PmzrE9EEDOGtgAOlgfBo_AK0poaXY/edit#gid=238579369)
      - Need to add oldest mice to collection list (including non-REL mice)
      - This Friday (13th March), collecting aged TRAP and all single-housed REL brains for OCT (15 mice) *16th: Done, in -20*
    - LCM:
      - Write github protocol for freezing, sectioning, and staining that allows brightfield visualisation of TH +ve cells
      - Call Carl Zeiss to confirm delivery time. PO not received - asked finance for update *16th: PO received, waiting on delivery estimate*
      - Sectioning with Greg Daubney this Thursday: 10-12:00:
        - I have 4 frozen brains: 
          1:  TRAP114.1c OVX unknown, cerebellum removed and placed against base of MB
          2:  TRAP111.1f OVX +ve, cerebellum removed and placed against base of MB
          3:  REL126.1g OVX -ve, half cerebellum cut, rostral brain cut, placed against base of STR, to section cerebellum into MB first
          4:  TRAP111.3f OVX -ve, half cerebellum cut, rostral brain cut, placed against base of STR, to section cerebellum into MB first
    - Spatial transcriptomics:
      - Write github protocol for freezing, sectioning and optimisation experiment from ST manual
      - Write github protocol for preparing oligo-dt slides, as described in Surmodics codelink manual
    - Long-read sequencing:
      - Do nanopore cDNA RT of OVX+ and OVX- TRAP IP MB RNA and use UTR primers to amplify full length transcript.
        - Using RNA from MB IP of:
          1:  REL112.3d OVX +ve
          2:  REL113.3a OVX -ve
        - Used 9 ul RNA in RT reaction. Performed 25 cycle, 2 min extension PCR with PRM primers and LongAmp Taq.
        - Will then quantify and use x ng cDNA for UTR PCR
        - Need to aliquot IDT VNP and SSP
      - Primers to internal exons have been ordered:
        - Perform 'RACE' amplification, with touchdown PCR when these have arrived (Thursday)
        
  
## 6th March 2020
- This week, performed two TRAPs:
  - First TRAP (3rd March) using 4 mice, new LYS/HS buffers with bought in KCl stock. Switched to P300 multichannel for HS washes and performed 6x 300 ul washes in 600 ul total HS volume. Problem of extra nonspecific RNA appears to be solved.
  - Second TRAP (5th March) using newly prepared bead matrix: Followed Heiman protocol, using 1 hour PL incubation and 1.5 hour antibody incubation. Made 4x matrix in 1x 1.5 ml protein lo-bind tube. Made new LS buffer with Sigma NP-40 for this. Used 1X matrix for experiment: 0.125x for MB, 0.0625x for DS/VS. Striatal RNA reliably detected across all samples (see ribogreen).
  - So plan to make more matrix and test function in 1/2 mice before applying to cohort.
- Meeting with NCR/GH:
  - Meeting with RWM requested for week of 16th
  
## 2nd March 2020
- Week before last, performed TRAP on 4x REL mice, 0.125x MB, 0.0625x DS/VS, 6 dounces, 12 samples, 96 well plate wash, p1200 multichannel, new buffers made up in a batch to cover all TRAP collections, had to use SIGMA NP-40 and made up KCl 2M and sterile filtered before using. Ribogreen looked positive for successful isolation of TRAP RNA and washing of beads, with more in VS than DS for all samples. 
- Last week, performed batches 1 and 2 of sample collection. At least two problems:
  - Some samples contained far too much RNA. This could be due to inadequate washing, or a property of the new bead matrix? This always happened in triplicate, also (ie. MB, DS and VS of one animal contained far too much RNA), so could this also be due to the lysis buffer? Fluorochem NP-40 is used in the lysis buffer, not SIGMA.
  - Some samples contained no/almost no RNA in the DS/VS samples. This could be a problem with washing/bead matrix/buffers.
- Plan is to switch to using a P20-P300 12-channel pipette, using 600 uL HS buffer per wash with the pipette set at 300 ul. I will switch to performing 6 washes, as per 2018, and will thoroughly mix before magneting. I will also make up new LS/LYS/HS buffers. On Tuesday 3rd, I will do TRAP with 4 mice, MB/DS/VS, 2 mice using 0.125x MB, 0.0625x DS/VS, 2 mice on 0.25x MB, 0.125x DS/VS with the P300 multichannel and 6 washes to evaluate.
- The issue of no/very little RNA requires comparison with a new bead matrix. If Tuesday/Wednesday's experiment still shows samples with no RNA, I will make 3x of new bead matrix (using new LS buffer) for comparison.
  
##
- So far have collected MB/DS/VS ELISA samples from 8 REL and 1 TRAP mouse. Have ordered new Abcam GFP ELISA kit.
- Performed TRAP on 3 REL mice, MB/DS/VS with matrix series: (MB: 0.5, 0.25, 0.125x, DS: 0.25, 0.125, 0.0625x, VS: 0.25, 0.125, 0.0625x).
  - Obtained RNA of good RIN from all samples and performed Smartseq2 with 25 cycles on 2.3 ul of elute (not normalised before smart seq). Yielded ~ 2-500 ng cDNA for qPCR after
  - qPCR shows greater enrichment in TH/DAT than 2018 and potentially some depletion of GAD/GFAP (calculations to follow)
  - 96 well plate works well, but needs optimising in weighing down so it doesn't float in water
  
- Wednesday 19th plan to trial a 1 mouse per condition TRAP. 6 dounces, rinse between. Whichever matrix concentration appeared optimal (likely smallest possible).
  
## 11th February 2020
- Have requested 1 young REL, 1 old REL and 1 old TRAP mouse for ELISA with original striatal dissection and comparison of cortex at SNpc zone and frontal zone. Will perform incubation overnight as written for 7th January 2020.
- Setting up for TRAP, have taken inventory of major items. Have stock for 24x matrix.
- Preparing buffers:
  ```
  (Volumes for 50 ml)
  Low-salt buffer:
    - 20 mM HEPES KOH (1 ml 1 M)
    - 150 mM KCl (3.75 ml 2 M)
    - 10 mM MgCl2 (500 ul 1 M)
    - 1% NP-40 (1x 10 % ampoule)
    - H2O (39.75 ml)
    + 1 Protease inhibitor tablet per 10 ml at use
    + 0.5 mM DTT at use
    + 100 ug/ml CHX at use
    + 10 ul/ml rRNasin and Superasin at use
    
   Tissue-lysis buffer:
    - 20 mM HEPES KOH (1 ml 1 M)
    - 150 mM KCl (3.75 ml 2 M)
    - 10 mM MgCl2 (500 ul 1 M)
    - H2O (44.75 ml)
    + 1 Protease inhibitor tablet per 10 ml at use
    + 0.5 mM DTT at use
    + 100 ug/ml CHX at use
    + 10 ul/ml rRNasin and Superasin at use
   
   High-salt buffer:
    - 20 mM HEPES KOH (1 ml 1 M)
    - 350 mM KCl (8.75 ml 2 M)
    - 10 mM MgCl2 (500 ul 1 M)
    - 1% NP-40 (1x 10 % ampoule)
    - H2O (34.75 ml)
    + 0.5 mM DTT at use
    + 100 ug/ml CHX at use
  
   Dissection buffer:
    - 1x HBSS (5 ml 10x)
    - 2.5 mM HEPES KOH (125 ul 1 M)
    - 35 mM Glucose (1.75 ml 1 M)
    - 4mM NaHCO3 (200 ul 1 M)
    - H2O (42.93 ml)
    + 100 ug/ml CHX at use
    
  Per TRAP, per MB, DS, VS x1, require:
    - 3 ml Dissection Buffer
    - 6 ml Low-salt Buffer for matrix
    - 4 ml Lysis Buffer
    - 20 ml High-salt Buffer
  
  So use 100 ul rRNAsin/Superasin per.
  ```
  - Preparing 3x matrix:
  ```
  1.  Thaw 3x of each mAb
  2.  Resuspend dynabeads by vortexing for > 30 sec
  3.  Wash 3x in PBS
  4.  Wash in solution A (0.1 M NaOH, 0.05 M NaCl, in DEPC-H2O)
  5.  Wash in in solution B (0.1 M NaCl in DEPC-H2O)
  6.  Resuspend biotinylated protein L in 500 ul PBS (to 1 ug/ul) and aliquot into 4x 120 ul volumes and freeze at -20 *C
  7.  Incubate 9 mg (900 ul original volume) dynabeads with 360 ul biotinylated protein L for 1 hour with end over end rotation
  8.  Perform 5x washes with PBS 3 % BSA sterile-filtered
  9.  Incubate beads with 3x mAb room temperature overnight, end over end rotation
  10. Wash 6x in LS buffer, 1 ml per wash
  11. Add 0.03% proclin 300 to the final product
  12. Store in an Eppendorf Protein LoBind tube at 4*C
  ```

## 7th February 2020
- This week have piloted using Sch1 hood. Took MB, DS, VS, CTX (from adjacent to midbrain) from 3 young, 3 old mice (OVX+/-) to obtain samples for ELISA. Ran ELISA this evening after 3 hour incubation. Signal in MB, VS and CTX, but not in DS. May require longer incubation (overnight) and also different CTX region for control sample. Alternatively, take sample from TRAP (non-REL) animal. Standard curve also did not work well. Do in duplicate next time and take care with washes to use fresh tips when washing standards? Performed 5x 350 uL washes today. Looks like wash buffer is PBS-Tween.
- Will re-run ELISA next week using Nunc plates. 
  - Before re-running, take MB, DS, VS, CTX from 1 young, 1 old REL + 1 young TRAP mouse. Take CTX from very exterior frontal region, before dissecting other tissues.
- Need to order more TMB (check small fridge next to antibody freezer for TMB)
- SWH has ordered 1 ml Eppendorf DNA loBind 96-well plates. Need to test these with TRAP sample handling.
- SWH has ordered a manual 8x 1 ml multichannel pipette. 
- So next major stage is to decide on appropriate DS/VS dissection (to obtain signal from both striatal samples)
  - Then perform trap on MB with 1/2x, 1/4x, 1/8x matrices, STR with 1/4x, 1/8x, 1/16x from 3 mice
- Dounce definitely required for TRAP: Pestle does not sufficiently homogenise tissue
  - So need 12x dounces 

## 30th January 2020
- Taking DS/VS (new approach) and SNpc/VTA as described on 29/01/2020 from REL119.4e. Showing Guusje Haver how to extract RNA in parallel. Aim to extract, cDNA synth with VILO IV and qPCR today. Add new SNpc/VTA primers tomorrow when they arrive.
- In parallel, further optimisation of R code to generate plots for OPDC report. At this stage, would like to show WGCNA module information from 18m MB OVX and 3' UTR lengthening in STR.

## 29th January 2020
- Run qPCR of first batch of dissected MB/STR.
  - DS/VS dissection improved. 5-fold enrichment of Coch in DS and 1.5-fold in Ace, compared to 2.3-fold Coch and 1.2-fold in former dissection method. SNpc/VTA not succesfully with Calb1. Repeat using dissection picture shown:
  ![](https://github.com/peterkilfeather/pk_trap/blob/master/dissection/SNPC_VTA_dissection_2.jpg)
    
    Description as in [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3487291/):
    > **Dissection Technique for isolating the midbrain dopaminergic neuropil. Technique for dissecting substantia nigra (SN) from ventral tegmental area (VTA).** 1. Remove overlying cortex and hippocampus (these structures have already been removed in the above image). 2. Make a vertical cut separating the pigmented area of the SN from the VTA. 3. Make a horizontal cut at or near the midline just above the dorsal and lateral most part of the SN. 4. Tease the SN away from the rest of the brainstem. Once SN has been removed, similarly tease the VTA away from the rest of the midbrain. Note: this section is typically observed at coordinates, relative to Bregma, at AP 5.7.

- Used Allen Brain Atlas to identify additional genes with either SNpc or VTA enrichment, relative to surrounding structures. [Calb1](https://mouse.brain-map.org/experiment/show/79556672), [Kcns3](https://mouse.brain-map.org/experiment/show/77371817), [Adrbk2](https://mouse.brain-map.org/experiment/show/74724787) and [Sdc2](https://mouse.brain-map.org/experiment/show/72080155) are most different. Ordered primers:
  ```
  Adrbk2 Origene
    FW: CGGACAAACTCTGCTTCATCCTG
    RV: CGCTGGCATAAAACCGCATCTC
  Kcns3 Origene
    FW: TGCTTACCTGCCACTCTGAGGA
    RV: CAGTTCCTCCATCACATGCAGC
  Sdc2 Origene
    FW: GAACAGAGCTGACATCCGATAAG
    RV: GGGATGTTGTCAGAACTGGACTC
  ```

## 28th January 2020
- RNA extracted. See lab book for concentrations. 0.5 ug cDNA made, waiting for calb1 primers.
- Preparing OPDC report, so going through TRAP R project and optimising.
- Switching to use updated kallisto quantification of PK1 and KW2 samples (November alignment, with updated human chromosome 4).
- New pseudoalignment percentages:
  ```
  for i in */run_info.json ; do cat $i | grep p_pseudoaligned | cut -d ' ' -f 2 | cut -d ',' -f 1 >> p_pseudoaligned.txt ; done
  for i in */ ; do echo ${i%/} ; done > sample_codes_in_order
  paste sample_codes_in_order p_pseudoaligned.txt > pseudoalignment_stats.txt
  ```
- Optimising R Project to generate report on 2019 progress for OPDC report.
  
## 27th January 2020
- Extract RNA for SNpc/VTA, DS/VS comparison. Use [ONT extraction method](https://github.com/peterkilfeather/pk_trap/blob/master/Extracting%20human%20RNA%20with%20TRIzol.pdf). cDNA synth with VILO IV kit. 

## 24th January 2020
- Preparation for TRAP:
  - Considering dissection of SNpc/VTA and Dorsal/Ventral STR. Testing dissection of SNpc/VTA using Calb1 primers to compare, as in [Gao, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23638090).
    - Primers using [Origene design](https://www.origene.com/catalog/gene-expression/qpcr-primer-pairs/mp201535/calb1-mouse-qpcr-primer-pair-nm_009788)
    ```
    FW: CTTGCTGCTCTTTCGATGCCAG
    RV: GTTCCTCGGTTTCGATGAAGCC
    ```
    - Dissecting TRAP111.3e SNpc/VTA CP/NAcc for testing Calb1/Coch/Ace differences. This is replicate number 1. Tissue is stored in box label 'IBIS'. Will need to extract RNA, cDNA synth and run qPCR.
  - In parallel, will need tissue for eGFP ELISA. Once qPCR shows region enrichment, move to collect tissue for ELISA.
    

## 9th January 2020
Notes on [Caffrey, 2007](https://www.ncbi.nlm.nih.gov/pubmed/17555970) and [Beevers, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28689993):
  - Chromosome 17q21 contains 1.3-1.6 Mb in LD, including MAPT locus. The LD is due to a 900 kb chromosomal inversion thought to be due to non-allelic homologous recombination between long coding repeats flanking the region: Two haploytype families result: H1 and H2.
  - Alternative splicing of MAPT exons 2 and 3 determines the number of N-terminal inserts. Splicing of exon 10 determines the number of microtubule binding repeats (3R or 4R).
  - Although MAPT has as strong association with PD, there is no tau tangle pathology.
  - H2 haplotype expresses 2x as much exon 2+3+ transcript as H1. H1 expresses 40% more exon 10+(4R) transcript than H2, with regional differences, suggesting a link between regional vulnerability and polymorphism-driven splicing change.

## 20th November 2019
- Awk to convert sample metadata into yaml config for PAQR_KAPAC:
  ```bash
  awk 'BEGIN {FS=","} {if($11="b2") print $0}' /zfs/analysis/pk_trap/r_master/metadata/sample_metadata_110719.csv | awk   '{if($2>249195) print $0}' | awk '{if($2<281000) print $2 ": {bam: " $2 ", type: " $5 "}"}'
  ```
  
## 15th November 2019  
- [Zhao, 2019](https://www.mdpi.com/1422-0067/20/1/212/htm): **Take a look at this review**:
  - May be worth running experiments in midbrain to identify translation start sites, ribosome footprinting...
- Remapping KW_DAT samples using new index: ensembl mouse 98 cdna + ncrna + all ensembl homo_sapiens SNCA transcripts (as opposed to using the NCBI sequence, as was done previously). Also generated modified GTF with human SNCA added, annotated as originating from 'human4' chromosome. cDNA fasta updated with 'human4' name in `/stripe/references/mouse/kallisto/ens98`. Used `--genomebam` option to generate BAM alignments to genome: Produced modified FASTA with entire mouse primary assembly + human chromosome 4, labelled as 'human4', for viewing in IGV:
  ```bash
  while read s ; R="" ; for i in *$s* ; do R+="$i " ; done ; do kallisto quant -i /stripe/references/mouse/kallisto/ens98/ens98_mmus_cdna_ncrna_hSNCA.idx -o /zfs/analysis/kw_trap/KW_DAT-TRAP/kallisto/kallisto_ens98_hSNCA_20191115/$s -t 8 --bias $R -b 100 --genomebam --gtf /stripe/references/mouse/kallisto/ens98/Mus_musculus.GRCm38.hsap_98_20191115.gtf.gz ; done >> /zfs/analysis/kw_trap/KW_DAT-TRAP/kallisto/kallisto_ens98_hSNCA_20191115/log.txt 2>&1 < kw_batch_2_codes
  ```
- To split a multi-feature fasta into its constituent features:
  ```bash
  csplit -s -z /path/to/INPUT.FA '/>/' '{*}'
  for i in xx* ; do n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; mv "$i" "$n.fa"
  done
  ```
- Then plan to test Sleuth for estimating transcript abundance changes in KW_DAT dataset.
  
## 12th November 2019
  - Priority to identify 3' UTRs differentially expressed between axon and cell body, or age, or OVX.
  - Two routes to identifying 3' UTRs from our data: de novo reconstruction, or using annotations including PolyAsite and Gencode. Both routes only capture a low percentage (10-30%) of what is actually there, and there is low agreement between these methods ([Chen, 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz068/5522019)).
  - A third strategy would require long read sequencing to capture full-length transcripts, including the 3' UTR. This would likely require an enrichment for genes of interest, to capture low abundance isoforms. The enrichment could either be probe-based, for shorter transcripts, or long-range PCR for longer transcripts.
  - According to [PolyASite DB](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz918/5588346):
    - Alternative first exons (dictated by choice of promoter) and alternative terminal exons contribut most to variation between human transcript isoforms ([Reyes, Huber, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5778607/))
    - Proliferating cells express short 3' UTR transcripts due to 3' end processing at coding region-proximal poly(A) sites (Sandberg), and vice versa for resting, differentiated cells
    - RNA binding proteins bind to 3' UTRs, so the length of this region can influence expression, through stability, translation, localisation, etc.
    - Poly(A) site databases have enabled discovery of novel polyadenylation signals, novel isoforms, such as 'intronic' polyadenylation in immune cells [Singh, 2018](https://www.nature.com/articles/s41586-018-0465-8)
    - ML in condition-dependent poly(A) site usage led to identification of RNA binding protein modulators [Gruber, 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1415-3)
    
  - In [PolyASite DB](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz918/5588346), they avoid spurious internal priming of poly(T) primers and clustered sites that are closely spaced and share regulatory signals (because these are likely due to imprecision of the 3' end data). They provide overall and per-sample usage quantification of the poly(A) site in each cluster. They only used samples sequenced from total RNA. The full workflow is described [here](https://github.com/zavolanlab/polyAsite_workflow).
  - AJ Gruber and M Zavolan's [review](https://www.nature.com/articles/s41576-019-0145-z.pdf):
    - Important point reiterated: transcript can have the same coding sequence, but different UTR, leading to different effect.
    - At the 3' end of most RNA pol II derived transcripts, endonucleocytic cleavage occurs at a poly(A) site, defined by a specific motif. This is followed by addition of a poly(A) tail. Most human genes have multiple poly(A) sites and there is condition-dependent regulator expression/binding to sequence or structures in pre-mRNA, influencing which poly(A) site is used. Therefore alternative cleavage and polyadenylation (APA) can alter the CDS or 3' UTR.
    - 3' UTRs ~140 bp in worm to 1-2 kb in human. 
    - 3' processing incolves over 80 proteins, of which ~ 20 are core proteins:
      - The core consists of four subcomplexes: **Cleavage and polyadenylation factor (CPSF), Cleavage stimulation factor (CSTF), Cleavage factor I (CFI) and Cleavage factor II (CFII) + some other proteins, including Symplekin and poly(A) polymerase (PAP).**
    - AAUAAA is the canonical poly(A) signal. This PAS is present at most PA sites, ~21 nucleotides upstream of the cleavage site. The PAS is recognised by the CPSF.
    - There are other variants of AAUAAA functioning as PAS and they have a position-dependent profile, also.
    - Core PAS regulatory motifs are enriched at distal poly(A) sites compared to proximal sites for tandem poly(A) sites in protein-coding genes (tandem meaning the proximal PAS gives rise to a truncated form of the distal PAS 3' UTR). A variant of AAUAAA (AAUACA) does not show differential enrichment, reflecting the sequence-specific efficiency in guiding cleavage and polyadenylation.
    - 3' end cleavage occurs ~21 nucleotides downstream of the PAS, typically adjacent to an adenosine. CPSF3 performs endonucleocytic cleavage, and the PAP addas a poly(A) tail to the cleavage product. Poly(A) binding protein 1 (PABPN1) bind to the nascent poly(A) tail, preventing the interaction between BPSF and PAP when the tail reaches ~250 nucleotides (limiting poly(A) tail length). Poly(A) tail lengths do differ and influence translational efficiency.
    - QAPA has the drawback that low expression isoform quantification accuracy is low and the majority of experimentally determined poly(A) sites are not associated with annotated transcript isoforms (although I think QAPA tries to overcome that by merging novel sites with annotated isoforms...)
    - In differentiation, genes with a single 3' UTR undergo changes in mRNA abundance, whereas genes with multiple 3' UTRs undergo changes in isoform ratio to achieve tissue specificity.
    - Long and short 3' UTRs are over-represented in the soma and neurites of neurons, respectively [Ciolli, 2019](https://academic.oup.com/nar/article/47/5/2560/5258023)
    - Activation of monocytes, B cells and effector T cells involves 3' UTR shortening [Sandberg2008](https://www.ncbi.nlm.nih.gov/pubmed/18566288). The NFATC1 transcription factor accumulkates upon T cell activation, due to possessing a short 3' UTR, triggering effector T cell activation. The switch from membrane bound to secreted IgM in B cell activation is due to 3' end processing of IgM at an APA site downstream of an intermediate exon, leading to formation of mRNA with a composite terminal exon.
    - 25% of experimentally-determined poly(A) sites in the DB are located in introns, because the corresponding isoforms have not yet been annotated, although some recently have [Singh2018](https://www.nature.com/articles/s41467-018-04112-z)
    - As well as tandem UTRs being processed at proximal and distal poly(A) sites, isoforms can be generated by PAS located at intronic sites, changing the protein that is expressed by truncation of the coding sequence. Alternatively, a PAS could be located in an alternative terminal exon that can be subject to alternative splicing: When spliced in, the transcript is cleaved at this PAS, when spliced out, the canonical PAS at the typical terminal exon is used (cassette terminal exon (TE)).
    - Typically, in an individual cell, only one possible 3' UTR isoform is detected [Gruber2018kapac](https://www.nature.com/articles/s41592-018-0114-z.pdf). Therefore, correlating poly(A) patterns within individual cells and the expression level of mRNAs encoding regulatory proteins could lead to identifying additional factors determining poly(A) site selection.
    - miRNA-mediated repression is slightly biased towards long 3' UTR transcripts, although this difference is small, as miRNA repression is most pronounced at binding sites close to 3' UTR boundaries. It has been shown that 3' UTR shortening can potentiate repression where there are miRNA binding sites located in the center of the long 3' UTR isoform.
    - RBP binding sites are short and 'degenerate': Many proteins interact at the same site at different stages in the mRNA life cycle.
    - In disease, changes in sequences elements can lead to loss or gain of alternative polydenylation, change 3' end processing and change mRNA expression.
    
- [Gruber2018tectool](https://www.nature.com/articles/s41592-018-0114-z.pdf), predicted novel isoforms containing alternative terminal exons:
  - Used the Ensembl [transcript level support](https://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html) scheme to compare predicted sequences with the annotation to see if there was some evidence. Most predicted isoforms had Ensembl evidence. 
  - Predictions were longer than Stringtie and Cufflinks and contained a higher percentage of canonical PASs ~21 nucleotides upstream of the PAS. TecTool also had greater inter replicate agreement than these other tools. 
  - In a multi-tissue experiment, they found novel isoforms that were the most expressed transcript of the corresponding gene. 
  - In a 201-T Cell single cell dataset, they found that the abundance of predicted isoforms matched the range of annotated isoforms. Considering only reads that splice into the 5' splice site of the terminal exon, once the TPM reaches 1-2, the isoform is detected in multiple cells. However, multiple isoforms were rarely present in a cell at the same time. 
  - For a tissue bulk RNA seq dataset, they generated updated annotations for each sample, then merged the annotations of replicate samples of each tissue, to obtain a tissue-specific annotation. They do not describe how they merge. They then used Salmon to quantify transcripts.


    
    
    - The mechanisms governing 3' UTR length are not fully characterised. Knockdown of individual core 3' end processing factors can both shorten and lengthen 3' UTRs
    - The ratio of isoforms could be caused by selective degradation. Increased short 3' UTR abundance in spermatogenesis is due to selective degradation of one long UTR isoforms by Tudor domain-containing protein 6 and regulator of nonsense transcripts 2
    
    
    
    
    
  - There are different types of polyA site:
    > In order of decreasing priority: TE, terminal exon; EX, exonic; IN, intronic; DS, 1,000 nt downstream of an annotated terminal exon; AE, anti-sense to an exon; AI, anti-sense to an intron; AU, 1,000 nt upstream in anti-sense direction of a transcription start site; IG, intergenic
  - Snca in the human, has > 20 polyA sites described ([PolyASite](https://polyasite.unibas.ch/search)). 
  
## 6th November 2019  
  - REL breedings - started (see personal google drive for spreadsheet).
  - Redo PCR of REL115.3a-f, 115.2a, 117.3e, 117.2c, 123.1i, 122.1f, 122.1j
## 28th October 2019
- [Google Slides of Figures](https://docs.google.com/presentation/d/1DEyPslen4YshVAra9fRdYCM5koILknd5MapfkIGl8G0/edit?usp=sharing)
- Done three comparisons of DS/VS vs MB (2-fold, 4-fold and 8-fold cut-offs). The axon vs MB = all DS vs MB genes that are also VS vs MB, not DS+VS vs MB DESeq2 comparison. Only included genes that are IP vs Input enriched. At the 2-fold striatum vs MB level, 11/107 PD GWAS genes (from [Nalls, 2019](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwjl-cGB777lAhXyQkEAHYROAcAQFjAAegQIAhAB&url=https%3A%2F%2Fwww.biorxiv.org%2Fcontent%2F10.1101%2F388165v2&usg=AOvVaw3j3ES6ys5-3TjxHKVay13h)) are dorsally enriched vs MB. 8 are ventrally enriched vs MB. CHD9, CTSB, KPNA1 are the dorsal only enriched GWAS genes.

## 23rd October 2019
- Marina Chekulaeva: RNA localisation in neurons. protein transport, mrna localisaiton and local translation efficiency  (3 mechs to subcellular protein localisaiton). Microporous membrane to separate mESC neurons body from downward growing axons/dendrites.Neurite localised mRNA and proteins (range >1 to >2.5 log2 FC. 
- Zipcode: cis-element in RNA in 3' UTR to mediate localisation. Recruit RBP to tether to motor protein and regulate localisaiton or translation stability/efficiency. Also miRNA binding sites and RNA modifications
- Nzip protocol: neurite soma separation protocol: pool of tiled RNA oligos to tile 3' UTRS ; make cDNA, clone into vector library, transfect and RNAseq to identify which fragments mediate neurite RNA localisation.
- Out of 20000 detected isoforms, 4000 are neurite enriched
- Cdc42 GTPase: cytoskeletal dynamics: 2 isoforms in humans with alternative last exon (exon 6 or exon 7). Two proteins produced: Different functions; 7 has more axon specification, 6 more soma. Last ten amino acids differ. mCherry targeted E7 utr transfected to check localisation. E6 and more soma than neurite (reversed 7). Puro-PLA: visualise local translation: van Dieck 2015). 
- RNA affinity capture + SILAC to take the two 3' UTRs, incubate with lysates and measure lights peaks to identify proteins that preferentially bind each 3' UTR isoform. Qki, Ptbp2, two splicing factors have binding affinity to E7. 
- Localised RBPs: Proteomics data of neurites/soma: find 29 RBPs enriched in neurites. Some known to be important for RNA localisation, so others suspected to be. PAR-CLIP of localised RBP: sequence the RNA bound to the RBPs. Also generate RBP knockouts and see how RNA distributions changes in soma/neurite separation model. 
- Nxf7 knockout increases neurite expression of several synaptic genes, including Snca. Soma not as affected. 
- PAR-CLIP with flag tagged Nxf7 for identifying all targets.
- MiRNAs in RNA localisation and local translation stability/efficiency: miRNA expression analysis of neurites/soma. 150 that are more present in neurites than soma. AGO-CLIP (argonaut) to IP argonaut, sequence what is bound - use dominant negative form of protein that interacts with argonaut and recruits polyadneylation complex: use short peptide to bind argonaut but not bind ccr4 = dominant negative block of miRNA expression. Then increase expression of RNA hook to relieve block. AGO-hook works better than AGO antibody for CLIP. 
- Mutations in RBPs: SMN1, TDP-43, FUS: RBPs that are implicated in motor neuron disorders. Induce these mutations and identify factors that induce changes in RNA localisation across multiple models. JPND consortium
  
## 1st-10th October 2019
- Work on establishing pipeline to quantify 3' UTR:CDS ratio for protein coding genes with constitutive exons. Validated Sudmant results succesfully, showing D1 neuron-specific aging change in ratio. Application of pipeline to TRAP MB/STR data shows now major global change. Worth investigating differential ratios of genes between conditions
- To identify axon/dendrite localisation signals, have been reading [Comprehensive catalog of dendritically localized mRNA isoforms from sub-cellular sequencing of single mouse neurons](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-019-0630-z#MOESM7).
  - To quantify isoform-level expression:
    - Use the last 500 nt of 3' end of gene as isoform quantification feature. This normalises length differences between 3' UTRs ( + the majority of their reads mapped within 500 nt of the 3' end). 3' ends less than 500 nt apart were merged into a single quantification feature, resulting in non-overlapping features. Genes with only one expressed 3' isoform were removed from analysis. In their case, most reads fell within 500 nt of the 3' end, so for genes with a 3' UTR >500 ng in length, these criteria avoid counting reads mapping to regions where the coverage will be low because of 3' bias. In our TRAP data, the 3' bias is..., so an alternative strategy is to simply map to each 3' UTR and calculated the length-normalised coverage of said feature.
    - Count how many 3' UTRs are present in the annotation per gene
    - Count how many 3' UTRs are expressed (reads > 0 in x/total samples)
    - Calculate length-normalised coverage of each 3' UTR
    - Calculate total number of reads across all 3' UTRs per gene
    - Calculate fraction of total reads per 3' UTR
  
## 29th September 2019
- Worked on getting 'ARMOR' (Soneson) snakemake pipeline working for RNA QC/mapping. Definitely need to look at RNASEQQC package for alternative splicing saturation analysis.
- Command to rsync file to EC2 instance:
  ```bash
  rsync -Pav Mus_musculus.GRCm38.dna.primary_assembly.ERCC92.hSNCA_20190924.fa -e "ssh -i /home/peter/.ssh/nanopore.pem" ubuntu@ec2-18-130-156-211.eu-west-2.compute.amazonaws.com:/nanopore/
  ```

  
## 26th September 2019
- Work has focussed on building a python pipeline to count reads mapping in each GTF interval for CDS and three prime utrs of filtered genes (described in 20th September section).
- Currently have number of reads per sample per interval. For each sample, need to condense all CDS into one CDS per gene with a total length and total number of reads. Same goes for three prime utr. Any genes without a three prime utr will be excluded. Then the coverage for each region will be multiplied by (1/length) and the UTR value will be divided by the CDS value and this output will be summarised to log base 2.

## 20th September 2019 
- awk commands to select protein coding genes with one annotated stop codon from GTF file: 21718 protein coding genes to start, of which 7736 have one annotated stop codon:
```bash
awk -F "\t" '$3 == "stop_codon" { print $9 }' Mus_musculus.GRCm38.97.gtf | grep 'gene_biotype "protein_coding"' | tr -d ";\"" |  awk '{print $2}' | sort | awk -F " " '{gene_counter[$1] += 1} END {for (gene_name in gene_counter){print gene_name, gene_counter[gene_name]}}' | sort | wc -l
#21718
awk -F "\t" '$3 == "stop_codon" { print $9 }' Mus_musculus.GRCm38.97.gtf | grep 'gene_biotype "protein_coding"' | tr -d ";\"" |  awk '{print $2}' | sort | awk -F " " '{gene_counter[$1] += 1} END {for (gene_name in gene_counter){print gene_name, gene_counter[gene_name]}}' | sort | awk '$2 == 1 {print $1}' > single_stop_codon_genes.txt
```

- In terms of exons, if you select 'ensembl_havana' supported transcripts there are 217388 'ENSMUSE...' exons:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -F "\t" '{print $9}' | tr -d "\";" | grep -oE "ENSMUSE[0-9]*" | wc -l
  #217388
  ```
  same story if you awk for exon in 3rd column:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | awk '$3=="exon"' | wc -l
  #217388
  ```
  In the ENSMUSE annotation, 186266 are uniquely named, so some exons are listed more than once with the same ENSMUSE ID: (note a regexp has been used with awk's match to return the matching ENSMUSE[0-9]\*, because the exon ID field moves positions in the 9th column of GTF files)
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' match($0, /ENSMUSE[0-9]*/) { exonName=$1":"$4":"$5":"$7; print(exonName, "\t", substr($0, RSTART, RLENGTH))}' | sort | awk '{print $2}' | sort | uniq | wc -l
  #186266
  ```
  If you take exons by coordinates, there are 217388:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' match($0, /ENSMUSE[0-9]*/) { exonName=$1":"$4":"$5":"$7; print(exonName, "\t", substr($0, RSTART, RLENGTH))}' | sort | awk '{print $1}' | sort | wc -l
  #217388
  ```
  of which 185344 are unique:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' match($0, /ENSMUSE[0-9]*/) { exonName=$1":"$4":"$5":"$7; print(exonName, "\t", substr($0, RSTART, RLENGTH))}' | sort | awk '{print $1}' | sort | uniq | wc -l
  #185344
  ```
  So there are more unique ENSMUSE IDs than actual coordinate exons. 
  If you list all exons by coordinates and ENSMUSE IDs and find duplicate coordinate rows, there are 922, explaining the difference:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' match($0, /ENSMUSE[0-9]*/) { exonName=$1":"$4":"$5":"$7; print(exonName, "\t", substr($0, RSTART, RLENGTH))}' | sort | uniq | awk '{print $1}' | sort | awk 'c[$1]++; c[$1]>=2' | sort | tail
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' match($0, /ENSMUSE[0-9]*/) { exonName=$1":"$4":"$5":"$7; print(exonName, "\t", substr($0, RSTART, RLENGTH))}' | sort | uniq | awk '{print $0}' | grep X:73351809:73351923:+
  #X:73351809:73351923:+    ENSMUSE00001211455
  #X:73351809:73351923:+    ENSMUSE00001214437
  ```
  Build a list of constitutive exon IDs:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; match($0, /ENSMUSE[0-9]*/); printf("%s\t%s\t%s\t%s\n",exonName,geneName,transcriptName, substr($0, RSTART, RLENGTH)); }' | sort | uniq | awk -v FS='\t' '{ eCount[$4]++; tCount[$3]++; exonHost[$4]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { print $4 }}' > constitutive_exons.gtf
  wc -l constitutive_exons.gtf
  #168908 constitutive_exons.gtf
  ```
  and non-constitutive exon IDs:
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'transcript_source "ensembl_havana"' | grep 'ENSMUSE' | awk -v FS='\t' ' { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; match($0, /ENSMUSE[0-9]*/); printf("%s\t%s\t%s\t%s\n",exonName,geneName,transcriptName, substr($0, RSTART, RLENGTH)); }' | sort | uniq | awk -v FS='\t' '{ eCount[$4]++; tCount[$3]++; exonHost[$4]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]!=gCount[exonHost[i]]) { print $4 }}' > non_constitutive_exons.gtf
  wc -l non_constitutive_exons.gtf
  #17358 non_constitutive_exons.gtf
  ```
  Build a GTF containing excluding multi-stop codon genes and non-constitutive exons: (Just filtering for single stop codon genes removes non-constitutive exons)
  ```bash
  cat Mus_musculus.GRCm38.97.gtf | grep 'gene_biotype "protein_coding"' | grep -f single_stop_codon_genes.txt  | grep -vf non_constitutive_exons.gtf > mus_ensembl_97_protcod_singlestop_constitexons.gtf
  ```
## 19th September 2019 
- Meeting with RWM, NCR. Points:
  - Spiking TRAP -ve lysate with GFP to determine GFP-dependent non-specific binding
  - Doing a qPCR for 3 synaptic genes from axonal TRAP
  - 3' UTR truncation in aged mice TRAP: How do they truncate? At a stop codon? Is aSyn truncated?
  - Get a list of calcium channels with interest for alternative splicing
  - Michael Coleman: Look up Nmnat2 and Sarm1 for Wallerian Degeneration relevance
  - **Meet again in 3 weeks (11AM 11th October) with update on:**
    - Oxidative truncation
    - Alternative splicing
    - Motif enrichment
    
## 18th September 2019 
- Set up REL breeding for EM perfusion: MB and STR

## 11th September 2019 
- Aim to calculate coverage before and after stop codon of transcripts in midbrain and striatal TRAP data. Compare young and old coverage, and see if it is further changed in OVX.

## 10th September 2019 
- [Sudmant, 2018](https://www.ncbi.nlm.nih.gov/pubmed/30485811) paper demonstrates age-related 3' UTR accumulation using TRAP and Gtex.
  - Steps to use:
    - Mapped using STAR followed by CIRCexplorer to quantify circular RNAs. STAR parameters:
      >-chimSegmentMin 15 -chimJunctionOverhangMin 15 parameters 
    - Generated a "termination codon ratio", defined as the log base 2 of the ratio of the mean read coverage after the annotated stop codon to the mean read coverage before the annotated stop codon.
    - Genes were filtered to use constitutive exons of non-overlapping genes, excluding genes with multiple annotated protein coding stop codons.
    - For plotting coverage upstream and downstream of the stop codon, coverage was normalized to windows of 1000 bins per gene (with gene structure placed below plot). The coverage value for the window with lowest coverage (per gene) was then subracted from every window and the resulting coverage was normalised to sum to 1, log transformed and smoothed using a Gaussian kernel 100 windows wide using the smth function in R smoother package. This could be done in tidyverse?
    
## 9th September 2019
- Meeting with RWM Friday. Agenda:
  >1. Biological follow up of current TRAP data. I need to be able, on Dec 3rd at the Wellcome Trust, to demonstrate in some way that at least one exemplar of the TRAP data in the Wellcome proposal is “true”. The TRAP data is novel and exciting, but I still think a bit vulnerable to miserable naysayers.
  >2. New TRAP cohort
  >3. Kathie’s paper.
  
- [Nalls, 2019](https://www.biorxiv.org/content/10.1101/388165v3) is the latest PD GWAS, identifying 90 signals. This list could be used for prioritising TRAP genes to pursue. GBA was included, but was in 'Hardy-Weinberg Disequilibrium'. 
  >Under Hardy-Weinberg assumptions, allele and genotype frequencies can be estimated from one generation to the next. Departure from this equilibrium can be indicative of potential genotyping errors, population stratification, or even actual association to the trait under study  [Turner, 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/). It has been consistently noted that many more SNPs are out of HWE at any given significance threshold than would be expected by chance. SNPs severely out of HWE should therefore not be eliminated from the analysis, but flagged for further analysis after the association analyses are performed.
  
- Have run 4 dissociation trials so far, with the most recent trial following the [Brewer, 2007](https://www.ncbi.nlm.nih.gov/pubmed/17545985) protocol. Obtained 3.2 million cells from one midbrain on last trial, with >75% viability. Cultured by overriden by bacteria so disposed of the plate. After talking to Nora and Siv, the likelihood of succesfully plating and surviving dopamine neurons is low, so switch back to flow analysis. Discussion with Natalie, recommended to use fixation and antibody staining as done for iPSCs, so stain for Th, in case the endogenous eGFP and FFN102 dye are insufficient. Charmaine is sending protocol and will plan experiment for Wednesday 11th September.
  - Initial plan: Use THTR surplus mice. 2-3+ mice, 5 weeks old. Brewer protocol, with D-AP5, Kynurenic Acid, EDTA for Ca2+-free, 5% Trehalose, 37&deg;C papain stage, ventilate HibA to oxygenate for (?) hours before protocol (discuss). Will use cortex tissue for Th-ve control region. 



## 30th July 2019
- Reversing order of notes (newest first)
-  Following lab meeting, priorities are:
  - Finish THTR paper:
    - qPCR for Cited2, Gfra1, Id2, Foxp2, Kitg showed only Foxp2 as significantly differentially expressed between ctrl and triplication. Emailed Jimena and she supplied the list of triplication DEGs that are Foxo3 targets ![](https://github.com/langkilfeather/pk_trap/blob/master/foxo3_trip_rnaseq.png). Ttr and Nnat are at the top of list when considering basemean of expression. Ordered primers based on Origene sequences for Ttr, Nnat, Slc39a7 and Cadm1. Depending on Ct values from previous qPCRs, dilute cDNA 1/2 or 1/4 to have enough for testing.
  - Get direct RNA seq going:
    - See minion repo
  - Dissociation

## 8th July 2019
- Searched for cre and eGFP sequences in IP sample reads
- Used Trinity to generate transcript sequences de novo
- eGFP-L10a: 'TRINITY_DN144_c0_g1_i8' length 1849
- IRES-Cre: 'TRINITY_DN144_c0_g5_i1' length 1805
- Cleared zfs snapshots:
  ```
  for snapshot in `sudo zfs list -H -t snapshot | cut -f 1` ; do sudo zfs destroy $snapshot ; done
  ```
- Created 'sample_codes.txt': cat the files, uniq, cut sample code, pipe into file:
  ```
  for i in fastq/*.gz ; do echo $i | cut -d '_' -f 3 >> temp ; done
  cat temp | sort | uniq >> sample_codes.txt
  ```
- Kallisto quant to quantify transcripts: 
  ```
  while read s ; do echo "$s" ; R="" ; for i in *$s* ; do R+="$i " ; done ; echo kallisto quant -i ../trinity_out_dir/trinity.idx -o kallisto/$s -t 8 --bias $R ; done < sample_codes.txt
  ```
- Summarize percentage pseudoaligned: 
  ```
  for i in */run_info.json ; do cat $i | grep p_pseudoaligned | cut -d ' ' -f 2 | cut -d ',' -f 1 >> p_pseudoaligned.txt ; done
  ```
- Generate ensembl_ERCC mouse index, annotation provided at mapping: 
  ```
  STAR --runMode genomeGenerate --runThreadN 8 --genomeDir Mus_musculus.GRCm38.dna.primary_assembly.ERCC --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.ERCC.fa
  ```
- Run STAR with on the fly annotation and 2 pass mapping with no shared memory: 
  ```
  while read s ; do echo "$s" ; R1="" ; R2="" ; for i in fastq/*${s}_1.fastq.gz ; do R1+="$i," ; done ; for j in fastq/*${s}_2.fastq.gz ; do R2+="$j," ; done ; STAR --genomeDir /references/Mus_musculus.GRCm38.dna.primary_assembly.ERCC/ --runThreadN 8 --outFilterType BySJout --outFilterMultimapNmax 20 -- alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --genomeLoad NoSharedMemory --outFileNamePrefix star_pc_ercc/${s}_ --readFilesCommand zcat --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --sjdbGTFfile /mnt/blue/reference/Mus_musculus.GRCm38.97.gtf --sjdbOverhang 74 --genomeSAindexNbases 14 --readFilesIn $R1 $R2 ; done < sample_codes.txt
  ```
- Imported unbound mb/ds/vs counts and plotted deseq normalised values for eGFP-L10a and Cre. Cre not detected in unbound ds/vs
- Imported ip mb/ds/vs counts and plotted deseq normalised. Cre not detected in ip ds/vs
- Cutadapt trimming of PK_TRAP samples improves kallisto mapping by 0.5% with initial trinity transcriptome
- Cutadapt trimming of PK_TRAP samples improves STAR mapping by 1% with ensembl GRCm38 genome

## 9th July 2019

**Why is there enrichment of astrocyte markers:**

![](https://github.com/langkilfeather/pk_trap/raw/master/astro_dopa_axon_plot.png)
![](https://github.com/langkilfeather/pk_trap/raw/master/cre_ip_ub.png)
- *In axon samples:*
  - Slc6a3 is 80-fold enriched, basemean 541
  - Th is 12-fold enriched, basemean 4231
  - Gfap is 1.93-fold enriched, basemean 6883
  - Aldh1l1 is 3-fold enriched, basemean 6436
  - S100b is 4-fold depleted, basemean 1014
  - eGFP-L10a is not enriched
  - Cre is enriched, but too low to quantify (input samples are 0, ip is 6-40)
- *In midbrain samples:*
  - Slc6a3 is 24-fold enriched, basemean 56462
  - Th is 66-fold enriched, basemean 247150
  - Gfap is 2.7-fold depleted, basemean 5803
  - Aldh1l1 is 5.5-fold depleted, basemean 767
  - S100b is 6-fold depleted, basemean 1681
  - eGFP-L10a is 4-fold enriched
  - Cre is 15-fold enriched

- Slc6a3 and Th are enriched in axonal samples, confirming specific binding of dopamine marker transcripts. Cre is enriched, and qPCR needs to be done for accurate measure.

- However many astrocytic markers are also enriched. This is either due to 
  1. ectopic Cre expression **or** 
  2. non-specific binding of striatal transcripts
  3. exchange of ribosomes/cDNA between neighbouring cells (Dougherty/Sakers, 2017)
  4. coprecipitation of astrocyte ribomes/RNA (Dougherty/Sakers, 2017)

- A TRAP experiment using 4 Cre-ve mice yielded no RNA. In contrast, the yield of RNA from Cre+ve striatum ranges between 15 and 40 ng (which is close to the range obtained from midbrain (50 ng - 100 ng)).

- In specifically bound, cell-type specific RNA, I assume there should be clusters of functionally related genes. Conversely in non-specifically-bound RNA, I assume that genes should have a weaker functional relationship to one another.

- Downloading and preparing Sakers, 2017 astrocyte TRAP data, to compare IP/Input of 3 bio replicates [paper](https://www.pnas.org/content/114/19/E3830) 
  - EBI file parsing and aspera download:
  ```
  while read s ; do awk '{ print $11 } ' >> links.txt ; done < PRJNA300571.txt
  while read s ; do awk '{ print $10 } ' >> md5sum.txt ; done < PRJNA300571.txt
  awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' accessions.txt | awk -F ";" 'OFS="\n" {print $1, $2}' | awk NF | awk 'NR > 0, OFS="\n" {print "ascp -QT -l 10m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh" " " $1 " ."}' > aspera.txt
  cat aspera.txt | parallel " {} "
  ```
  - Fastqc looks odd: GC bias in some samples and variable read duplication rate
  - Using GRCm38.ERCC92.cdna for quantification with kallisto as done for my samples
  - Kallisto: `for i in *.gz ; do kallisto quant -i /references/Mus_musculus.GRCm38.cdna.ncrna.ERCC92.idx -o kallisto/${i%.fastq.gz} -t 8 --bias --single -l 200 -s 30 $i ; done` [Using -l 200 -s 30 for single end](https://www.biostars.org/p/252823/). 78% mapping so acceptable.
  - Added DESeq2 chunk

- Interestingly, Sakers paper uses this method to determine enriched transcripts:
  >*PAP-enriched transcripts were determined by the intersection of transcripts enriched in SN input/cortex input [false-discovery rate (FDR) < 0.1, log2 (fold-change) > 1.5] and transcripts enriched in PAP TRAP/SN input [FDR < 0.1, log2 (fold-change) > 1.5]. These transcripts are listed in Dataset S2.*

  >*PAP-depleted transcripts were determined by the intersection of transcripts depleted in SN input/cortex input [FDR < 0.1, log2 (fold-change) < 1.5] and transcripts enriched in cortex TRAP/cortex input [FDR < 0.1, log2 (fold-change) > 1.5]. These transcripts are listed in Dataset S3.*

- Cre and eGFP-L10a primers ordered:
  > - IRES-Cre_TRAP_qPCR_PK_fw: ACCTGTTTTGCCGGGTCAGA
  > - IRES-Cre_TRAP_qPCR_PK_rv: TCCAGGGCGCGAGTTGATAG
  > - eGFP-L10a_TRAP_qPCR_PK_fw: CTGTATTAGCCCGGGCCCTC
  > - eGFP-L10a_TRAP_qPCR_PK_rv: GTATGGGTACATGGTGGCGTC

- Produced a table in master_08.07.19 to evaluate expression of striatal cell type markers in axon_ip_vs_ub results table. Astrocyte/OPC/microglia rank top.

## 10th July 2019

- Creating merged fastq files for trinity transcriptome generation. Taken 6 fastq files, 2 MB, 2 DS, 2 VS (3/18 month):
  ```
  294132
  295120
  296108
  215118
  216106
  217189
  while read s ; do echo $s ; for i in trimmed/*${s}*1.fq.gz ; do echo -n "$i " >> R1.txt ; done ; for j in trimmed/*${s}*2.fq.gz ; do echo -n "$j " >> R2.txt ; done ; done < trinity_input_codes
  zcat $(cat R1.txt) >> r1.fq
  zcat $(cat R2.txt) >> r2.fq
  pigz *.fq
  Trinity --seqType fq --max_memory 16G --left r1.fq.gz --right r2.fq.gz --CPU 8 --output trinity_10.07.19
  ```
- Also get medium spiny neuron dataset? PCA of highest enriched genes to compare axon samples with astro and msn?

## 11th July 2019
- Rebuilding WGCNA for midbrain. Need to assess suitability of batch 3 samples for joining with batch 2 (KW).
  - Need to decide gene filtering threshold for input
- Trinity construction had to be run on CGAT cluster. Tried reducing in silico normalisation to 30 from 200 on hudson but still not enough

## 12th July 2019
- Testing Cufflinks for de novo transcript assembly
  ![The pipeline on the right is the current recommended cufflinks workflow](http://cole-trapnell-lab.github.io/cufflinks/images/tuxedo_workflow.png)
  - Input bam files should be sorted. Testing using the 6 bam files sorted on 10th July:
  ```
  mkdir sorted
  while read s ; do echo $s ; for i in ${s}_Aligned.out.bam ; do samtools sort -@ 8 -o sorted/${i}_Aligned_sorted.out.bam $i ; done ; done < ../trinity_input_codes
  ```
## 13th July 2019
  - Cre and eGFP-L10a primers arrived. Will run primers on the cDNA used for checking DAT/TH enrichment. Testing primers on RNA extracted from midbrain tissue:
## 14th July 2019
  - **Meeting with RWM, NCR**
    - THTR: Finish final qPCRs, manuscripts, prepare cover letter:
      - Triplication RNA quantified with ribogreen. All < 1ng/ul. Going to run SMART-SEQ2 with 2.3ul of RNA regardless of concentration (not diluting samples before amplification).
      - Why it is ‘fab’
      - PLOS, not so much about novel in cover letter, but thorough, in depth, good science
      - But it is the first TRAP profile of PD genetic model
      - Check AD/Huntingtons TRAP papers? Otherwise this is first genetic neurodegeneration TRAP paper
  - Filter microglia/macrophage/msn by Friday
  - Report genes of interest status for grant by Friday:
    - Synapto-Everything
    - Dynamin 1,2,3
    - Micu1
    - Calcium channels
    - Kif1a/b
    - Sortin/Nexin1
  - Find examples of S100b -ve astrocytes to explain depletion of S100B in axon-TRAP samples
  - Find out how many cell soma enriched transcripts are found in Sakers Aldh1l1 TRAP dataset (particularly Slc6a3)
  - For future, controls for ectopic Cre expression should be:
    - Take 1 whole cortex from GFP+ve mouse and run TRAP
    - Repeat 4x striata from GFP-ve mice
  - Produce confident Venn diagram of axon/mb by Friday
  - Caleb wants to know axonally translating TFs: Make list
  - **Next meeting Friday afternoon. Time TBC**

## 18th/19th July 2019
- Not a clear enrichment of MSN markers:
  ![](https://github.com/langkilfeather/pk_trap/blob/master/msn_markers.png)
- Enrichment of cholinergic interneuron markers:
  ![](https://github.com/langkilfeather/pk_trap/blob/master/cholinergic_markers.png)
- Enrichment of predominantly tissue macrophage markers. Important microglia markers are not enriched. Markers taken from [Barres paper]([Microglia and macrophages in brain homeostasis and disease, 2017](https://www.nature.com/articles/nri.2017.125#microglia-and-brain-macrophages))
  ![](https://github.com/langkilfeather/pk_trap/blob/master/microglia_macrophage_markers.png)
- Meningeal and choroid macrophages:
![](https://github.com/langkilfeather/pk_trap/blob/master/tissue_macrophage_markers.png)

- Using McKeever, 2017 [Cholinergic neuron gene expression differences captured by translational profiling in a mouse model of Alzheimer's disease.](https://www.ncbi.nlm.nih.gov/pubmed/28628896):
  ```
  awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' PRJNA387000.txt | cut -f9 | awk -F ";" 'OFS="\n" {print $1, $2}' | awk NF | awk 'NR > 1, OFS="\n" {print "ascp -QT -l 20m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh" " " $1 " ."}' > download.txt
  cat download.txt | parallel " {} "
  ```
  - Fastqc/Multiqc:
    ```
    mkdir qc
    fastqc -t 8 -o qc *
    multiqc -o qc qc/
    ```
  - Kallisto: `for i in *1.fastq.gz ; do kallisto quant -i /references/Mus_musculus.GRCm38.cdna.ncrna.ERCC92.idx -o kallisto/${i%.fastq.gz} -t 8 --bias $i ${i%1.fastq.gz}2.fastq.gz ; done`

- Ayata, 2018 to examine microglia involvement: [Epigenetic regulation of brain region-specific microglia clearance activity](https://www.nature.com/articles/s41593-018-0192-3):
  - Subsetted to remove single nuclei files
  -  Added ip/input and region columns
  ```
  awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' PRJNA387000.txt | cut -f9 | awk -F ";" 'OFS="\n" {print $1, $2}' | awk NF | awk 'NR > 1, OFS="\n" {print "ascp -QT -l 20m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh" " " $1 " ."}' > download.txt
  cat download.txt | parallel " {} "
  ```
  - Fastqc/Multiqc:
  ```
    mkdir qc
    fastqc -t 8 -o qc *
    multiqc -o qc qc/
    ```
  - Kallisto: `for i in *1.fastq.gz ; do kallisto quant -i /references/Mus_musculus.GRCm38.cdna.ncrna.ERCC92.idx -o kallisto/${i%.fastq.gz} -t 8 --bias $i --single -l 200 -s 30 $i ; done` Single-end: Alignment %:

- Bone marrow derived macrophage ribotag. [Jackson, 2018](https://www.nature.com/articles/s41586-018-0794-7?WT.feed_name=subjects_translation#data-availability)
- Interesting microglia paper, also using IgG control ribotag. Check this dataset for the amount of enriched transcripts in IgG control. [Haimon, 2018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114001)

## 20th/21st July 2019
- Pipelined EBI/aspera and started making snakemake file for fastqc/multiqc/kallisto
- Sakers astro data does not have enough depth to detect Slc6a3, to see if it was enriched in their dataset. The [Boisvert astrocyte aging dataset](https://www.cell.com/cell-reports/pdf/S2211-1247(17)31848-X.pdf) did different regions (visual cortex, motor cortex, spinal cord, hypothalamus, cerebellum) so testing the hypothalamus samples to check for Slc6a3 (they have 25-60 million reads per sample).
- There is Th enrichment in the McKeever cholinergic dataset, but no Slc6a3: Th expression is confirmed in GABAergic interneurons [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4465985/)
- 6000 genes are enriched in the cholinergic dataset: 2000/3800 axonal genes are in common. 2745/4000 MB enriched genes are in common with cholinergic.
- 433 genes are enriched in astros, cholinergic, and midbrain DA: Should they be returned to the axon-enriched list? I think so
- 940 genes common between macrophage ribotag and axon trap. 
- Send RWM/NCR email with different venn diagram versions. Filtering cholinergic genes has the issue of including neuronal genes that could be shared between all neuron types. ![Suggested venn diagram](https://github.com/langkilfeather/pk_trap/blob/master/axon_enriched_astro_macro_filter.png)

## 23rd July 2019
- Made changes to boxplot multi code for RWM. Made borders and text bolder
- Reading [Cuffdiff 2 paper](https://www.nature.com/articles/nbt.2450#ref28):
  - Count uncertainty: Up to 50% of reads map ambiguously to different transcripts because in higher eukaryotes, alternative isoforms share large amount of sequence and many genes have paralogs with similarity. Therefore transcript counts are estimates.
  - Count overdispersion: Although we use a Poisson model to estimate variability, the variability in count data between replicates is more than what we would expect in a Poisson distribution. This overdispersion increases with expression [(Anders, Huber, 2010)](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?holding=npg&cmd=Retrieve&db=PubMed&list_uids=20979621&dopt=Abstract) and Anders proposed using the negative binomial distribution to control for it.
  - Cuffdiff 2 estimate uncertainty by calculating the confidence that a fragment is correctly assigned to the transcript: if the transcript has more shared exons with other isoforms, there will be lower confidence. Similarly if there are only a few assigned fragments, there will be low confidence. This confidence is summarised in a 'beta distribution' and the overdispersion of counts is modelled with the negative binomial distribution, and it 'mixes them together', resulting in a beta negative binomial distribution.
- Kruskal-Wallis tests that a continuous variable has the same mean across multiple groups [Peter Langfelder blog](https://peterlangfelder.com/2018/11/25/working-with-categorical-variables/)
- Need to get WGCNA code simpler and with **clear gene names**


 

