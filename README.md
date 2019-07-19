DPAG fly TRAP 
Kevin Talbot 'Shatra'
Jakob Scaber 
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










- Lab meeting discuss:
  - DESeq2 normalisation vs TPM/RPKM/FPKM
  - Trinity/Cufflinks de novo reconstruction to isoform detection
  - Homer motif analysis
  - WGCNA midbrain
  - Nanopore
  - FFN102/Dissociation status



[Data Managment plan](https://github.com/sr320/LabDocs/wiki/Data-Management)
