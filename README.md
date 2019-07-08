Notes

## 8th July 2019
- Searched for cre and eGFP sequences in IP sample reads
- Used Trinity to generate transcript sequences de novo
- eGFP-L10a: 'TRINITY_DN144_c0_g1_i8' length 1849
- IRES-Cre: 'TRINITY_DN144_c0_g5_i1' length 1805
- Created 'sample_codes.txt': cat the files, uniq, cut sample code, pipe into file
- Kallisto quant to quantify transcripts: `while read s ; do echo "$s" ; R="" ; for i in *$s* ; do R+="$i " ; done ; echo kallisto quant -i ../trinity_out_dir/trinity.idx -o kallisto/$s -t 8 --bias $R ; done < sample_codes.txt`
- Imported unbound mb/ds/vs counts and plotted deseq normalised values (log2 transformed) for eGFP-L10a and Cre. Cre not detected in unbound ds/vs
- Imported ip mb/ds/vs counts and plotted deseq normalised









[Data Managment plan](https://github.com/sr320/LabDocs/wiki/Data-Management)
