# denorm
Denormalization and DEG detection for GSE23561 dataset.

### Data Resources

- Files under `GPR/` are available in [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23561), named as `GSE23561_RAW.tar`.

- Files under `CSV/` can be generated using `samet.sh`.

- Raw version of `GPL.csv` file can downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10775).

- Raw version of `SERIES.csv` file is available in [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23561), named as `Series Matrix File(s)`.

- Raw version of `MAP.csv` file is [Human ReadyArray](http://www.microarrays.com/docs/Human_MI_ReadyArray_genelist.xlsx). It is also available in [Microarrays Website](http://www.microarrays.com/resources.php).

- Download `9606.protein.links.v10.5.txt.gz` file from [STRING database](https://string-db.org/cgi/download.pl?sessionId=dgQBX5PHNs5I&species_text=Homo+sapiens). Export to `IMPORT/` directory and name as `PLINKS.tsv`.