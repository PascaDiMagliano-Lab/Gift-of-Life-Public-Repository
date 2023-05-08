Donor scRNAseq raw sequencing files can be acquired from [here](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003229.v1.p1), and tumor scRNAseq raw sequencing files can be acquired from [here](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002071.v1.p1). 

Samples were aligned to hg38 using cellranger v6 on the University of Michigan HPC. scRNAseq samples count matrices can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229413), and spatial transcriptomics raw data can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226829)

## scRNAseq Samples
You can download the files using this [link](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229413/suppl/GSE229413_RAW.tar), move them to "data" folder and extract the .tar file into a folder with the same name (GSE229413_RAW)

or run the following command from terminal while being in the 'data' folder

```
cd </PROJECT/TO/DATA/PATH/IN/PROJECT/FOLDER>
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229413/suppl/GSE229413_RAW.tar --output GSE229413_RAW.tar
mkdir GSE229413_RAW
tar -xzvf GSE229413_RAW.tar --directory GSE229413_RAW
```
## Nanostring GeoMx Samples

Similarly, you can download the `.dcc` files using this [link](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226829/suppl/GSE226829_RAW.tar), move them to "data" folder and extract the .tar file into a folder with the same name (GSE226829_RAW). We'll also need to download the following files

* [pkc file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226829/suppl/GSE226829_Hs_R_NGS_WTA_v1.0.pkc.gz)
* [segment annotation file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226829/suppl/GSE226829_PreQC_segments_annotation.xlsx)

or run the following command from terminal while being in the 'data' folder

```
cd </PROJECT/TO/DATA/PATH/IN/PROJECT/FOLDER>
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226829/suppl/GSE226829_RAW.tar --output GSE226829_RAW.tar
mkdir GSE226829_RAW
tar -xzvf GSE226829_RAW.tar --directory GSE226829_RAW
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226829/suppl/GSE226829_Hs_R_NGS_WTA_v1.0.pkc.gz --output GSE226829_Hs_R_NGS_WTA_v1.0.pkc.gz
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226829/suppl/GSE226829_PreQC_segments_annotation.xlsx --output GSE226829_PreQC_segments_annotation.xlsx
```


