```bash

```
# Task1_Epigenomics
This is the first task to be deliverd into the Epigenomics course, as a part of the Masters in Omics Data Analysis imparted by Universitat de Vic - UCC. For this exercise, we have cloned and followed the informations available at Beatrice Borsari's Github repository called Epigenomics_uvic (https://github.com/bborsari/epigenomics_uvic/wiki). 

```bash
git clone https://github.com/bborsari/epigenomics_uvic
```

Before starting,  we will use Docker and run the container using this command:

```bash
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

Once inside the container, we can start to operate.

## Exercise 4
### Task 1:
Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

To complete this task, we will create three principal folders using mkdir: one for data, one for annotation and one for analyses. Inside each folder, more folders will be created in order to have all newly data well organised.

```bash
cd ATAC-seq
mkdir data
mkdir data/tsv.files
mkdir data/bed.files
mkdir data/bigBed.files
mkdir data/bigWig.files
mkdir analyses
mkdir analyses/peaks.analyses
mkdir annotation
```

### Task 2
Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.

To start, data was obtained from the same source as in the hands on explanation, but the filter used will be the one specified on the ennunciate. To download it, the following command will be used:

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&biosample_ontology.term_name=stomach&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&type=Experiment"
```

Once our data has succesfully been downloaded, we will briefly explore it's content using
```bash
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```
The information needed is located in column 1 (file accession), for tissue column 11 (Biosample_term_name) and column 23 (Experiment_target). We will first parse the metadata file to download the files in the corresponding folders bigBed.files or bigWig.files.
This command sequence filters lines from the file metadata.tsv based on specific criteria, processes and formats the output using awk, sorts it, and then writes the result to a file named bigBed.peaks.ids.txt in the analyses directory.
```bash
grep -F "bigBed_narrowPeak" metadata.tsv | grep -F "pseudoreplicated_peaks" | grep -F "GRCh38" | awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' | sort -k2,2 -k1,1r | sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
cat analyses/bigBed.peaks.ids.txt
ENCFF287UHP     sigmoid_colon
ENCFF762IFP     stomach

cut -f1 analyses/bigBed.peaks.ids.txt | while read filename; do wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; done
ls data/bigBed.files
ENCFF287UHP.bigBed  ENCFF762IFP.bigBed  md5sum.txt
```
When downloading files, it's always advisable to check their integrity by verifying their MD5 hash, a sort of digital fingerprint of the file. MD5 hashes can be computed with the command md5sum, so this command performs an integrity check on a set of files. It compares the original MD5 sums stored in metadata.tsv with the computed MD5 sums of the corresponding files and reports any discrepancies.

```bash
file_type="bigBed"; ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt; cat data/"$file_type".files/md5sum.txt | while read filename original_md5sum; do md5sum data/"$file_type".files/"$filename"."$file_type" | awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}'; done > tmp; mv tmp data/"$file_type".files/md5sum.txt; awk '$2 != $3' data/"$file_type".files/md5sum.txt
```
Since no output was obtained we can conclude that the are no differences between the files, they have a correct integrity and are verified.

### Task 3
For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions). Hint: have a look at what we did here and here.

We will start by dowloading the gencode.v24.primary_assembly.annotation file and unncompress the gtf.gz file:
```bash
wget -P annotation "https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz"
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz
```

Now it's time to convert the gtf annotation file to a BED format. Specifically, we will retrieve gene body coordinates of protein-coding genes (chr, start, end, strand), remove mitochondrial genes (i.e. those located on chrM) and move from a 1-based to a 0-based coordinate system.

```bash
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf | grep -F "protein_coding" | cut -d ";" -f1 | awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' | sed 's/\"//g' | awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed
```
The next step is to convert bigBed files format into bed files, using bigBedToBed as follows:
```bash
cut -f1 analyses/bigBed.peaks.ids.txt | while read filename; do bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed; done
```
Now, we were provided the list of promoters ([-2 kb, +2 Kb] from TSS) of protein-coding genes and we stored this file inside the annotation folder with the name /gencode.v24.protein.coding.non.redundant.TSS.bed. Following this line of code we'll conduct intersection analyses for both tissues. We'll extract the first two columns from the data, where the first column represents the filename and the second column represents the tissue. These values will be used to calculate the intersection and the results will be stored in a text file with the tissue name.
To obtain the number of peaks the command wc -l will be used too.
```bash
 cut -f-2 analyses/bigBed.peaks.ids.txt | while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u | sort -u > analyses/peaks.analyses/ATACpeaks.within.promoters."$tissue".txt; done
root@97de1c55167d:/home/micalball/epigenomics/epigenomics_uvic/ATAC-seq# wc -l analyses/peaks.analyses/*.txt
  47871 analyses/peaks.analyses/ATACpeaks.within.promoters.sigmoid_colon.txt
  44749 analyses/peaks.analyses/ATACpeaks.within.promoters.stomach.txt
  92620 total
```





```bash

```
```bash

```


```bash

```
```bash

```



