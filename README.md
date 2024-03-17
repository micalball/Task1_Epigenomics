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
To obtain the number of peaks intersecting promoter regions the command wc -l will be used.
```bash
 cut -f-2 analyses/bigBed.peaks.ids.txt | while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u | sort -u > analyses/peaks.analyses/ATACpeaks.within.promoters."$tissue".txt; done
root@97de1c55167d:/home/micalball/epigenomics/epigenomics_uvic/ATAC-seq# wc -l analyses/peaks.analyses/*.txt
  47871 analyses/peaks.analyses/ATACpeaks.within.promoters.sigmoid_colon.txt
  44749 analyses/peaks.analyses/ATACpeaks.within.promoters.stomach.txt
  92620 total
```
As can be seen above, we found a total of 92620 peaks inside promoter regions, with 47871 corresponding to sigmoid colon and 44749 to stomach.
Now, we are going to repeat the same step but selecting the peaks in the BED file that do not intersect with the protein-coding gene body regions and we will save the results into separate files for each tissue.

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt | while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed -v | sort -u > analyses/peaks.analysis/ATACpeaks.outside.gene."$tissue".bed; done
root@97de1c55167d:/home/micalball/epigenomics/epigenomics_uvic/ATAC-seq# wc -l analyses/peaks.analyses/*.bed
  37035 analyses/peaks.analyses/ATACpeaks.outside.gene.sigmoid_colon.bed
  34537 analyses/peaks.analyses/ATACpeaks.outside.gene.stomach.bed
  71572 total
```
In this case, we found 71572 peaks outside gene coordinates with 34537 for stomach and 37035 peaks for sigmoid colon tissue.



## TASK 5

### Task 1
Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results.

```bash
cd ../epigenomics_uvic
mkdir regulatory_elements
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
Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

For H3K27ac
```bash
grep -F H3K27ac ../ATAC-seq/metadata.tsv|grep -F "bigBed_narrowPeak"|grep -F "pseudoreplicated_peaks" |grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |sort -k2,2 -k1,1r |sort -k2,2 -u > analyses/H3K27ac.bigBed.peaks.ids.txt && cut -f1 analyses/H3K27ac.bigBed.peaks.ids.txt |while read filename; do wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; done
```

For H3K4me1
```bash
grep -F H3K4me1 ../ATAC-seq/metadata.tsv|grep -F "bigBed_narrowPeak"|grep -F "pseudoreplicated_peaks" |grep -F "GRCh38" |awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |sort -k2,2 -k1,1r |sort -k2,2 -u > analyses/H3K4me1.bigBed.peaks.ids.txt && cut -f1 analyses/H3K4me1.bigBed.peaks.ids.txt |while read filename; do wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"; done
```

Now we will perform the check as in Task 4

```bash
for file_type in bigBed; do ../bin/selectRows.sh <(cut -f1 analyses/*"$file_type".peaks.ids.txt) ../ChIP-seq/metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt && cat data/"$file_type".files/md5sum.txt | while read filename original_md5sum; do md5sum data/"$file_type".files/"$filename"."$file_type" | awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' ; done > tmp && mv tmp data/"$file_type".files/md5sum.txt && awk '$2!=$3' data/"$file_type".files/md5sum.txt ; done
```
No results received so there are no differences. On the other hand, the list of candidate distal regulatory elements is the following:
```bash
cat data/bigBed.files/md5sum.txt
ENCFF977LBD     be29636550527e36c4755ea036531e75        be29636550527e36c4755ea036531e75
ENCFF844XRN     de679228721fb4055aa1f657c77c21a6        de679228721fb4055aa1f657c77c21a6
ENCFF724ZOF     c87fefbf41de3d291fa1d340a26627f5        c87fefbf41de3d291fa1d340a26627f5
ENCFF872UHN     2207b7b3378df7776e7ecdc2aa1a5de0        2207b7b3378df7776e7ecdc2aa1a5de0
```
No it's time to convert our bigBed files to Bed, using the bigBedToBed function. We will apply it to H3K27ac and then H3K4me1.

```bash
cut -f1 analyses/H3K27ac.bigBed.peaks.ids.txt|\
while read filename; do
        bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

```bash
cut -f1 analyses/H3K4me1.bigBed.peaks.ids.txt|\
while read filename; do
        bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

The next step is to find intersects of H3K4me1 and H3K27ac where distal reguylatory regions overlap. We will use bedtools intersect function for that.

```bash
root@8e1b947ec41c:/home/micalball/epigenomics/epigenomics_uvic/regulatory_elements# cut -f-2 analyses/H3K4me1.bigBed.peaks.ids.txt | while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b ../ATAC-seq/analyses/peaks.analyses/ATACpeaks.outside.gene."$tissue".bed -u > analyses/peaks.analyses/"$tissue".H3K27ac.outside.bed; done
```
```bash
cut -f-2 analyses/H3K27ac.bigBed.peaks.ids.txt | while read filename tissue; do bedtools intersect -a data/bed.files/"$filename".bed -b ../ATAC-seq/analyses/peaks.analyses/ATACpeaks.outside.gene."$tissue".bed -u > analyses/peaks.analyses/"$tissue".H3K27ac.outside.bed; done
```
Now check if our files have been succesfully created:
```bash
ls analyses/peaks.analyses
sigmoid_colon.H3K27ac.outside.bed  sigmoid_colon.H3K4me1.outside.bed  stomach.H3K27ac.outside.bed  stomach.H3K4me1.outside.bed
```
An find the intersect for each tissue:
```bash
for tissue in stomach sigmoid_colon; do cut -f-2 analyses/H3K27ac.bigBed.peaks.ids.txt | while read filename; do bedtools intersect -a analyses/peaks.analyses/"$tissue".H3K27ac.outside.bed -b analyses/peaks.analyses/"$tissue".H3K4me1.outside.bed -u > analyses/peaks.analyses/overlap."$tissue".bed; done; done
```
Lastly, we count the peaks inside the files created.

```bash
wc -l analyses/peaks.analyses/overlap*.bed
  7367 analyses/peaks.analyses/overlap.sigmoid_colon.bed
  4342 analyses/peaks.analyses/overlap.stomach.bed
 11709 total
```

### Task 3
Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file  regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.
We will use the two files created in the last step of the previous task. A filter for the desired columns containing the crucial information is used and results are stored in a new file called chr1.regulatory_elements.starts.tsv

```bash
for tissue in sigmoid_colon stomach; do
  awk 'BEGIN{FS=OFS="\t"} $1=="chr1" {print $4,$2}' analyses/peaks.analyses/overlap."$tissue".bed > analyses/peaks.analyses/chr1.regulatory_elements.starts.tsv
done
```
Now we count:
```bash
wc -l analyses/peaks.analyses/*.tsv
528 analyses/peaks.analyses/chr1.regulatory_elements.starts.tsv
```

### Task 4
Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3').
Using the provided line, we created the file.

```bash
awk 'BEGIN{FS=OFS="\t"} $1=="chr1" {if ($6=="+"){start=$2} else {start=$3}; print $4, start}' ../ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed > analyses/gene.starts.tsv
```
Now use wc -l to obtain the number of genes.
```bash
wc -l analyses/gene.starts.tsv
2047 analyses/gene.starts.tsv
```

### Task 5
Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

A new file called get.distance.py in epigenomics_uvic is created with the provided script:
```bash
nano ../bin/get.distance.py
```

```bash
root@8e1b947ec41c:/home/micalball/epigenomics/epigenomics_uvic/regulatory_elements# python ../bin/get.distance.py --input analyses/gene.starts.tsv --start 980000
ENSG00000187642.9       982093  2093
```

We get the expected result!


### Task 6
For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:

```bash
for tissue in stomach sigmoid_colon; do   cat analyses/peaks.analyses/chr1.regulatory_elements.starts.tsv |   while read element start; do     python ../bin/get.distance.py --input analyses/gene.starts.tsv --start "$start";   done > analyses/peaks.analyses/"$tissue".gene_distances.tsv; done
```

```bash
head -n 1 analyses/peaks.analyses/stomach.gene_distances.tsv
ENSG00000237330.2       1074307 6841
```

### Task 7
Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.

```bash

```

