# Repeat Browser Data Processing

---

In this repo, we provide a pipeline to process the transposable elements(TE),
which is used for getting alignment statistics of transposable elements and 
save the processed data into a proper file format that suitable for our repeat browser. 

This pipeline includes a modified [iteres](https://epigenome.wustl.edu/iteres/) pipeline and 
a python script to convert the analysis into [zarr](https://zarr.dev/) format. 
The output zarr files can be uploaded into our [repeat browser](https://repeatbrowser.org/) for the visualization.


## Prerequisites

---

### 1. Compile the iteres
```bash
$ cd iteres
# then run a make command here
$ make
```
if the 

### 2. Install the required packages of python



