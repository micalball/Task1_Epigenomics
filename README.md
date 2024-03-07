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




```bash

```



```bash

```


