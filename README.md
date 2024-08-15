# DGAT-cancer
DGAT-cancer was a method for predicting cancer driver genes. It integrates the predicted pathogenicity of the somatic mutations (together with germline variants in the healthy population) with a topological network of gene expression in tumor tissues, and considering the expression levels in tumor and paracancerous tissues while employing 24 distinct features to predict novel cancer drivers.
## Algorithm
### 1.	Using Laplacian Score to select features
We used an unsupervised method, the Laplacian Score, to select features that have high power to preserve the local geometric structure of the features space. The detailed steps for the application of Laplacian Score in feature selection were as follows:  
(1) A network ![](https://latex.codecogs.com/svg.image?\mathit{N}) was constructed to connect all m candidate genes (![](https://latex.codecogs.com/svg.image?\mathit{N}\in&space;\mathbb{R}^{\mathit{m}\times&space;\mathit{m}})). For each pair of genes, we calculated the Euclidean distance, ![](https://latex.codecogs.com/svg.image?\left|\mathit{x}_{i}-\mathit{x}_{j}&space;\right|) between their feature vectors (![](https://latex.codecogs.com/svg.image?\mathit{x}_{i}) is the feature vector of gene ![](https://latex.codecogs.com/svg.image?\mathit{i}), and ![](https://latex.codecogs.com/svg.image?\mathit{x}_{j}) is the feature vector of gene ![](https://latex.codecogs.com/svg.image?\mathit{j})). Based on the Euclidean distance, if gene ![](https://latex.codecogs.com/svg.image?\mathit{i}) is among the top ![](https://latex.codecogs.com/svg.image?\mathit{k}) (here we set ![](https://latex.codecogs.com/svg.image?k=\left&space;[&space;0.01\times&space;m&space;\right&space;])) of the nearest genes to gene ![](https://latex.codecogs.com/svg.image?j), or the gene ![](https://latex.codecogs.com/svg.image?j) is among the top ![](https://latex.codecogs.com/svg.image?k) of the nearest gene to the gene ![](https://latex.codecogs.com/svg.image?i) (![](https://latex.codecogs.com/svg.image?i\neq&space;j)), we set the connection between ![](https://latex.codecogs.com/svg.image?i) and ![](https://latex.codecogs.com/svg.image?j) as ![](https://latex.codecogs.com/svg.image?N_{ij}=1). Otherwise, ![](https://latex.codecogs.com/svg.image?N_{ij}=0).  
(2) In the network ![](https://latex.codecogs.com/svg.image?N), if ![](https://latex.codecogs.com/svg.image?N_{ij}=1), we weight the edge by ![](https://latex.codecogs.com/svg.image?W_{ij}=e^{-\left\||&space;x_{i}-x_{j}\right\||^{2}) (![](https://latex.codecogs.com/svg.image?\mathit{W}\in&space;\mathbb{R}^{\mathit{m}\times&space;\mathit{m}})). Otherwise, ![](https://latex.codecogs.com/svg.image?W_{ij}=0). The weighted network reflects the local structure of the ![](https://latex.codecogs.com/svg.image?m) genes in the feature space.  
(3) Computing the Laplacian Score of each feature. Let ![](https://latex.codecogs.com/svg.image?y_{l}=\left&space;[&space;y_{l1},&space;y_{l2},&space;...,&space;y_{lm}&space;\right&space;]^{T}) denote the ![](https://latex.codecogs.com/svg.image?l)-th feature values for all m genes. We redefine ![](https://latex.codecogs.com/svg.image?y_{l}) by removing the mean from the samples as in Equation (1). The Laplacian Score (![](https://latex.codecogs.com/svg.image?LS)) is computed by Equation (2). According to the scores of each feature, we finally selected the top 20 features with the smallest ![](https://latex.codecogs.com/svg.image?LS) to be included in the prediction model for each type of cancer.  
![](https://latex.codecogs.com/svg.image?\widetilde{y}_{l}=y_{l}-\frac{y_{l}^{T}DI}{I^{T}DI}I,&space;\mathit{(1)})  
where ![](https://latex.codecogs.com/svg.image?D=diag(\sum_{j=1}^{m}W_{1j},\sum_{j=1}^{m}W_{2j},...,\sum_{j=1}^{m}W_{mj}),&space;I=\left&space;[&space;1,1,...,1&space;\right&space;]^{T}).  
![](https://latex.codecogs.com/svg.image?LS_{l}=\frac{\widetilde{y}_{l}^{T}L\widetilde{y}_{l}}{\widetilde{y}_{l}^{T}D\widetilde{y}_{l}},&space;\mathit{(2)})  
where ![](https://latex.codecogs.com/svg.image?L=D-W).  
### 2.	Data transformation  
In order to integrate multiple features of each gene into a risk score (Fig. 1), we combined the features by Hotelling and Box-Cox transformations, which converted the feature values into p-values. For a given scaled matrix of m genes with n features (P∈R^(m×n), here n=20), the Hotelling transformation is performed as: 
█(P^'=U(P-IM)^T,)
where U=[v_1;v_2;…;v_(n^' ) ]^T  with n^'≤n as the number of chosen principal components of P and V=[v_1;v_2;…;v_n ] are eigenvectors for the covariance matrix of P corresponding to decreasing eigenvalues with λ_1≥λ_2≥⋯≥λ_n. M=[1/m ∑_(i=1)^m▒〖P_i1,〗  1/m ∑_(i=1)^m▒〖P_i2,…,〗  1/m ∑_(i=1)^m▒P_in ]. Thus, transformed P^'∈R^(n^'×m). Then, the Box-Cox transformation is performed as follows, 
█(〖pi〗^'={█((〖pi〗^(β_i )-1)/β_i ,β_i≠0@log(pi),β_i=0)┤  ),
where pi is the i-th row vector of P^' and all elements of pi vector are forced to be positive before being transformed, and β_i is the parameter for transforming pi to 〖pi〗^'. 
Finally, we standardized each 〖pi〗^'  and calculated P-values for elements of 〖pi〗^'   that is in a standard Gaussian distribution. The P-values of the elements were combined as S(j) by Fisher's method (Fisher's combined probability test). S(j) is termed the score of gene j. 
S(j)=-ln⁡(〖_n'^2〗^(-1) (-2∑_(i=1)^(n^')▒log⁡(〖pi〗_j^' ) ))

### 3.	Gibbs sampling  
The scores (![](https://latex.codecogs.com/svg.image?S)) of genes were used as conditional probabilities in Gibbs sampling to obtain a convergent probability distribution of candidate genes. In the first round of sampling, Gibbs sampling was initiated by randomly selecting ![](https://latex.codecogs.com/svg.image?m') (![](https://latex.codecogs.com/svg.image?m'\leq&space;m)) genes from the candidate genes. The ![](https://latex.codecogs.com/svg.image?m') genes were assumed to have equal probabilities of being selected. Then, in the second round of sampling, another set of ![](https://latex.codecogs.com/svg.image?m') genes was sampled from the remaining ![](https://latex.codecogs.com/svg.image?m-m') genes weighted by their scores (![](https://latex.codecogs.com/svg.image?S)). The following rounds were the same as the second round. In each round, the selected frequency (the number of times a gene was selected divided by the total number of sampling) of each gene was updated. All the selected frequencies of candidate genes in the ![](https://latex.codecogs.com/svg.image?i)-th round were denoted as a vector (![](https://latex.codecogs.com/svg.image?m\times&space;1)), ![](https://latex.codecogs.com/svg.image?Freq_{i}). When the Euclidean norm of ![](https://latex.codecogs.com/svg.image?Freq_{i}-Freq_{i-1}) was smaller than ![](https://latex.codecogs.com/svg.image?E_{Gibbs}) (![](https://latex.codecogs.com/svg.image?E_{Gibbs}) was set as 0.01), the iteration was stopped. ![](https://latex.codecogs.com/svg.image?Freq_{last}) was assigned as the posterior probabilities (PP) of candidate genes.  
### 4.	Random permutation  
Then, we constructed a null distribution of PP in order to obtain the likelihood of a given gene being a cancer driver gene. The null distribution was generated by giving genes randomly weighted scores that were derived from the uniform distribution with the same range as the true scores. We generated 1,000,000 sets of null distributions of PPs for the m genes by running Gibbs sampling. For each gene, we obtained 1,000,000 random PPs. By counting the number of times that a random PP of a gene was larger than the real PP of the gene, an experience p-value was estimated. The p-values were adjusted by Bonferroni correction and we selected those genes with ![](https://latex.codecogs.com/svg.image?p_{adj}<0.01) as being significant.  
## Code running instructions  
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
