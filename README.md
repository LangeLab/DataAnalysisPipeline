# [Name of App here]
Here we illustrate the common workflow of [Name of app].
## Data Input
Datasets for input should be from the same set of samples, i.e., with the same meta data, and the first column of meta data should be unique IDs of samples. Upon datasets of intensities (protein, termini, peptide and other PTMs), the column of unique IDs for each row should be selected. These unique IDs will be used for annotation and data integration.

[Insert an example of datasets  here]

If any dataset of termini/PTM is available, the annotation tool will be available as well. Note that the information is retrieved from Uniprot service and [Name of app] only supports organisms  of human now. 

[Insert an example of annotation tool here]

If the dataset of peptide is available, the protein intensity calculator function will be available. The user can choose to sum up intensities of peptides in different ways to obtain protein intensities. Also the result can be used in following analysis.

[Insert an example of peptide calculator here]

## Quality Control
The standard QC focuses on the batch effect detection and correction on standard samples, which should be uploaded in this section. 
The sample QC section displays 5 plots of the datasets uploaded in Data Input section. If there is a column indicating replica, averaging across replica can be conducted at this step. Also if batch effect is detected in standard QC, correction can be conducted on the regular samples with settings inheritted from the standard QC.

[Insert an example of batch effect correction here]

## Pre-processing
This section is consist of filtering, normalization and imputation. Different settings can be conducted on different datasets and they are processed individually. 

- Filtering. Rows with NA's more than given percentage will be dropped.

- Normalization. Aside from median normalization, the user can choose normalization within sample grouping/protein accession for PTMs. 

- Imputation. "Down-shifted Normal Sampling" imputes the NA with random samples from a narrower Normal distribution with smaller mean, and calculation of such shifting and shrinkage is based on the sample grouping/condition the user choose. Other methods available are direct implementation of *impute* function of *MSnBase* package. 

## Statistical Inference
 This section mainly conducts statistical testing on differential expression of quantitiative intensity from each of all available datasets between different sample groups. Multiple testing p-value adjustment by B-H method is default. Further interactive settings include:

- Standardization.  Rescales the values into a range of [0,1]

- Blocking. Blocking is a technique for dealing with nuisance factors which are factors that have some effect on the response but of no interest. Hence the influence they have on the response variable needs to be minimized or ruled out. Currently the blocking factor can be 2-level factor only, 

-  Weighting. We introduced a function that allows the user to put different weights on original/imputed values in testing of differential expression. By specifying a rather small weight on imputed values may reduce the influence from imputation methods on statistical inference.

[Insert example here]

## Dimensionality Reduction

Dimensionality reduction techniques can help with problems such as large covariance matrix (which could bring trouble to computation), pairwise correlations between the variables to consider, or irrelevant variables prior to modelling. Valid methods including PCA, t-SNE and UMAP.

[Insert example here]

## Clustering
This section features unsupervised learning, inlucding hierarchical clustering, non-hierarchical clustering (k-means) and fuzzy clustering. Note that rows with NA values are not handled in this section and they are directly dropped. 

## Individual Protein Visualization
Selections of dataset to display are available in the Boxplot tab. Aside from displaying rows that are significant in differential expression, users can specify the IDs of rows to be displayed and IDs should be separated strictly with comma only. 

- Note that the domain structure and annotation PTMs is only available when "protein data" is selected in the boxplot tab and PTM dataset is available. The PTM dataset should include columns indicating protein accession, modification position and modifitation type. 

- To generate circosplot, at least one of termini/PTM dataset should be available along with protein dataset, and statistical inferences on datasets to include should be conducted in Statistical Inference section. 
