```bash

```
# Task1_Epigenomics
This is the first task to be deliverd into the Epigenomics course, as a part of the Masters in Omics Data Analysis imparted by Universitat de Vic - UCC. 
In this exercise, we have used and followed the informations available at Beatrice Borsari's Github repository called Epigenomics_uvic (https://github.com/bborsari/epigenomics_uvic/wiki). 
Before starting,  we will use Docker and run the following container using this command:
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

To start, data was obtained from the same source as in the hands on explanation, but the filter used will be the one specified on the ennunciate.  To download it, the following command will be used:

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&biosample_ontology.term_name=stomach&assay_title=ATAC-seq&biosample_ontology.term_name=sigmoid+colon&type=Experiment"
```




```bash

```



```bash

```



```bash

```


