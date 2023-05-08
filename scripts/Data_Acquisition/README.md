Donor scRNAseq raw sequencing files can be acquired from [here](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003229.v1.p1), and tumor scRNAseq raw sequencing files can be acquired from [here](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002071.v1.p1). 

Samples were aligned to hg38 using cellranger v6 on the University of Michigan HPC. scRNAseq samples count matrices can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229413), and spatial transcriptomics raw data can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226829)

You can download the files using this [link](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229413/suppl/GSE229413_RAW.tar), move them to "data" folder and extract the .tar file into a folder with the same name (GSE229413_RAW)

or run the following command from terminal while being in the 'data' folder

```
cd </PROJECT/TO/DATA/PATH/IN/PROJECT/FOLDER>
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229413/suppl/GSE229413_RAW.tar --output GSE229413_RAW.tar
mkdir GSE229413_RAW
tar -xzvf GSE229413_RAW.tar --directory GSE229413_RAW
```

