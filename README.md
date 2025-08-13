# Data Processing Pipeline for Small RNA Analysis in HtRNAPan

## Introduction 
  HtRNAPan, an integrated platform for comprehensive landscape analysis and functional annotation of herbal tRNAs. The database features: 1) Over 514,119 tRNAs and over 33,417,735 tsRNAs from 315 medicinal plant species. 2) Annotation of over100 tRNA modification enzymes across all species. 3) tRNA sequences with predicted modification types, positions, probabilities. 4) RNA modification isozyme profiles. 5) tsRNA sequences and regulatory roles in cross-kingdom. Users can intuitively explore tRNA panorama data via a user-friendly interface.


**cite from:**  (文章名称)[www.biorxiv.org]



## Pipeline Workflow
![pipline](img/pipline.png)

### 1. Data Download and Cleaning 

Retrieve datasets from the following sources:

- [NCBI](https://www.ncbi.nlm.nih.gov/)

- [Modomics](https://iimcb.genesilico.pl/modomics/)


  For Modomics data, we downloaded the original HTML files and performed data scraping using the Python script in the `01.data_scrape` directory as follows.

  ```sh
  #protein 
  python 01.data_scrape/01.modifications/run.py -i 01.data_scrape/01.modifications/test/1.html -o ./1.json
  
  #modification 
  python 01.data_scrape/02.proteins/run.py -i 01.data_scrape/02.proteins/test/1.html -o ./1.json
  ```

  

### 2.  tRNA 2D and 3D Structure Prediction

​	Predict 2D and 3D Structures Using the [VfoldPipeline_alone](https://rna.physics.missouri.edu/vfoldPipeline/index.html).

### 3.  tRNA Modification Prediction

<200b>  

### 4.  Analysis of Modification-Related Isozymes

​	Use 

### 5.  tsRNA Mining and Target Gene Prediction

<200b>  
