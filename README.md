# StellarPath: patient classifier <img src="[https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/sp/Logo.png](https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/sp/Logo.png)" width=350 align="right" />

## Brief description

StellarPath is a deep learning patient classifier which leverages state-of-the-art bioinformatics techniques to select biologically meaningful molecular and pathway features, the pathway space to integrate different omics and increase the software interpretability, the PSN paradigm to enable analysis and stratification of patients based on their pathway activity, and an artificial neural network to account for the similarity information and classify.

## Citation
StellarPath: Patient Classifier Integrates Multi-Omics in Patient Similarity Networks for the Detection of Stable Markers and Signature Pathways. Authors: **Giudice Luca, Mohamed Ahmed and Malm Tarja**

## Features

1. Free of hyper-parameters (input data are the same of common pathway analysis tools)
2. Analyses raw counts from their normalization to integration in pathways without need of user’s pre-processing (at the end, the user will have access also to the normalized counts)
3. Integrates multi-omics information in pathways to ease interpretability and not deviate from biology 
4. Uses a novel measure to determine how much two patients are similar. It considers how much and in which direction two patients regulate a specific pathway. 
5. Creates pathway specific PSNs where each node represents a patient is study, an edge connects two patients based on how much they are similar in a pathway.
6. Select molecules and pathways if their PSNs satisfy a criterium (one class has patients who are more similar than the members of the opposite class and the two classes are not similar, i.e. one class is cohesive and the other one is sparse).
7. Feeds the PSNs to a graph convolutional network that learns how to recognize the patients based on their similarities and how to predict the class of unseen individuals.
8. PSNs can be visualized with an ad-hoc implementation that reduces the number of patients which are visible as nodes without losing their information.

## Results

1. molecular markers (genes, miRNAs, …) which are differentially expressed, differentially stable and predictive of a patient’s class. 
2. pathways that are enriched by markers, significantly deregulated, represented by biologically meaningful PSNs, and predictive of a patient’s class. 
3. significant PSNs showing that the pathway activity of the patients of one class is precisely regulated with highly (positive/negative) stable molecules while the members of the opposite class present more variable and different expression levels.
4. patient’s centrality score that reflects how much each patient in study is similar to the members of its class and dissimilar from the others in a pathway. This enables the patient stratification. 

## Install

StellarPath works with R but depends on a python module.
You can install the R package with:

```r
devtools::install_github(
  repo="LucaGiudice/StellarPath",
  ref="main",
  dependencies = "Depends",
  build_vignettes = FALSE
)
```

Once the R package is installed, you can install the Python module.
The python version used to develop the method has been the 3.8.12.

- For non Windows-based systems, you need to install python 3 and conda/miniconda <https://conda.io/miniconda.html>.
- For Windows-based system, you need to use Visual Studio, install python 3 and conda/miniconda.

Once conda is installed, our R package include two functions to install a python enviroment linked to the package:

``` r
#For Windows-based system, run R as administrator to make installation work properly
library("StellarPath")
SP_install() #check conda, if conda is missing installs automatically miniconda, then it installs the python enviroment
SP_initialize(save_profile = TRUE) #activate the conda enviroment
```

**In case of any installation issue:**

We suggest to use the docker: ***lgiudice/py3_r4_bioinf:latest***. It has an Ubuntu OS that is suited to install and work with StellarPath.

## Vignettes

There are the following vignettes:

1. [Introduction to StellarPath's workflow, data and results](https://htmlpreview.github.io/?https://github.com/LucaGiudice/supplementary-files/blob/main/sp/Vignette_1.html)
2. [Introduction to StellarPath's how to plot the data](https://htmlpreview.github.io/?https://github.com/LucaGiudice/supplementary-files/blob/main/sp/Vignette_2.html)

## Workflow

<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/sp/workflow.png" />

Frame A shows the workflow:
  1. The method classifies two patients’ classes and works with both dense and sparse (e.g., somatic mutation) omics.
  2. It finds the significantly deregulated molecules and their enriched pathways.
  3. It determines how much each pair of patients is similar by comparing the values of the molecules belonging to a specific pathway. It uses the pathway-specific similarities to build a network. Only the PSNs, which show that one class is cohesive while the opposite one is not, are kept.
  4.  StellarPath receives an unknown patient, adds it to a significant PSN and predicts its class with a graph convolutional network. The latter step is repeated with all PSNs and the final prediction is made by consensus/major voting.
  5. StellarPath provides multiple output data: B) a PSN represents the pathway’s similarities due to differentially expressed and variant molecules, C) the enriched pathways are prioritized with the properties of their PSNs, D) both known and unknown patients are stratified, and quality checked based on how much are similar to the others in the enriched pathways.

## Integration system

StellarPath integrates omics with a knowledge-based approach in pathways and takes charge of both the ome's processing and results.
This means that, before putting a new ome available in the system, we will test the results in wet-lab to be sure that its integration with other data is not creating false biological findings. 
The current omics and supported integration is the following one:
<p align="center">
	<img src="https://raw.githubusercontent.com/LucaGiudice/supplementary-files/main/sp/System.png" width=55% height=55%/>
</p>

We will add other omics in the close future.

## License

BSD-3-Clause @ Giudice Luca
