# DGAT-cancer
DGAT-cancer was a method for predicting cancer driver genes. It integrates the predicted pathogenicity of the somatic mutations (together with germline variants in the healthy population) with a topological network of gene expression in tumor tissues, and considering the expression levels in tumor and paracancerous tissues while employing 24 distinct features to predict novel cancer drivers.
## Algorithm
### 1.	Using Laplacian Score to select features
We applied an unsupervised method called the Laplacian Score to select features that have a high ability to preserve the local geometric structure of the feature space. The process involved the following steps:
#### Network Construction: We constructed a network to connect all candidate genes. For each pair of genes, we calculated the Euclidean distance between their feature vectors. If a gene was among the top k nearest genes to another gene, or vice versa, we established a connection between them in the network.
#### Edge Weighting: In this network, we assigned a weight to each connection based on the Euclidean distance between the genes, using an exponential function. This weighted network represented the local structure of the genes in the feature space.
#### Laplacian Score Calculation: For each feature, we calculated its Laplacian Score, which measures how well the feature preserves the local structure of the data. Features with lower Laplacian Scores are considered more important. We selected the top 20 features with the lowest scores for inclusion in the prediction model for each cancer type.
### 2.	Data transformation  

To integrate multiple features of each gene into a risk score, we applied a combination of Hotelling and Box-Cox transformations, which converted the feature values into p-values. For a scaled matrix of m genes with n features (where n=20 in this study), we first performed the Hotelling transformation. This involved selecting a certain number of principal components from the feature matrix, based on the eigenvectors corresponding to the largest eigenvalues of the covariance matrix. Next, we applied the Box-Cox transformation to ensure that all feature values were positive, adjusting them based on specific parameters to enhance normality. Finally, the transformed features were standardized, and p-values were calculated, representing how each feature fits within a standard Gaussian distribution. These p-values were then combined using Fisher's method to produce a risk score for each gene.

### 3.	Gibbs sampling  
Gibbs sampling offered a refined probabilistic approach to identify and prioritize cancer driver genes by evaluating their posterior probabilities.

The scores (S) of genes were used as conditional probabilities in Gibbs sampling to obtain a convergent probability distribution of candidate genes.  Initially, a random subset of genes was selected from the candidate pool, with each gene assumed to have an equal chance of being chosen.  In subsequent rounds, another subset of genes was sampled from the remaining pool, but this time weighted by their scores.  This process was repeated, with the frequency of each gene being selected updated in each round.  The iteration continued until the difference between the selection frequencies in successive rounds was minimal.  The final selection frequencies were then assigned as the posterior probabilities of the candidate genes.

### 4.	Random permutation  
Then, we constructed a null distribution of PP in order to obtain the likelihood of a given gene being a cancer driver gene. The null distribution was generated by giving genes randomly weighted scores that were derived from the uniform distribution with the same range as the true scores. We generated 1,000,000 sets of null distributions of PPs for the m genes by running Gibbs sampling. For each gene, we obtained 1,000,000 random PPs. By counting the number of times that a random PP of a gene was larger than the real PP of the gene, an experience p-value was estimated. The p-values were adjusted by Bonferroni correction and we selected those genes with ![](https://latex.codecogs.com/svg.image?p_{adj}<0.01) as being significant.  

# Code running instructions  
1.	Organizing the feature file as “cancer.list.Rdata”. Gene symbols are in row and features are in column.  

2.	laplacian_and_hotelling_box-box_transform.R.   
The data file was then put into the laplacian feature selection to find the most important features. Then the selected features would be transformed using Hotelling and Box-cox to generate weights of genes, which were used as sampled weights in Gibbs sampling.  

3.	gibbs_sampling.R  
The weighted score file generated in Step 2 was used as input to obtain a convergent probability distribution of candidate genes. The final selected frequencies were assigned as the posterior probabilities (PP) of candidate genes.  

4.	random_permutation.R  
The weighted score file generated in Step 2 was used as input to construct a null distribution of PP to obtain the likelihood of a given gene being a cancer driver gene. By counting the number of times that a random PP of a gene was larger than the real PP of the gene, an experience p-value was given.   

5.	The p-values were adjusted by Bonferroni correction and those genes with ![](https://latex.codecogs.com/svg.image?p_{adj}<0.01) were selected as being significant.  
# Acknowledgement  
We thank for all the available mutation data and RNA-seq data for supporting this work. Somatic mutation data were downloaded from the Broad Institute GDAC Firehose Portal (http://gdac.broadinstitute.org/), the International Cancer Genome Consortium (ICGC) Data Portal (https://dcc.icgc.org/releases/current/Projects/) and The Cancer Genome Atlas (TCGA) (https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga). The germline mutation data of healthy individuals were collected from Phase3 of the 1000 Genomes Project (https://www.internationalgenome.org/data, GRCh38). RNA-seq data (RSEM normalized count, log2 transformed) of tumors from 12 cancer types in TCGA were downloaded from the UCSC Xena platform63 (http://xena.ucsc.edu/). The corresponding RNA-seq data of paracancerous tissues were downloaded from TCGA.  
# Note  
The pathogenic status of somatic and germline mutations were scored by 19 distinct approaches collected in dbNSFP (version: dbnsfp30a) of ANNOVAR, where each representing a different detrimental effect of mutations on gene function. Then we calculated the difference between the two profiles (annotation functional scores in tumor and healthy cohorts, respectively) by means of a unidirectional Earth Mover’s Difference (uEMD) score. These uEMD scores together with expression levels of genes in tumors and paracancerous tissues, scores generated by a topological gene expression network in cancer tissues.  
# Language
R (100%)
