# Log2 Fold Change and T-test Pipeline for Novel Molecular Targets Identification in Multiple Sclerosis 

## Overview

* MS is a chronic autoimmune disease that affects over 1 million people in the US alone, and 2.8 million people worldwide.
* Current MS therapies fail to prevent the progression of the disease, only slowing its progression
* Novel Research suggests that lipid peroxidation, as well as oxidative stress, may be possible drivers of axonal degeneration within lesions
* The objective of this research was to validate novel molecular targets for Multiple Sclerosis using statistical analysis and bio-omics tools
* This project is a custom-made Python pipeline that independently attempted verification of three molecular targets for MS - calprotectin, including S100A8 and S100A9, as well as ALOX15B, across 10 raw GEO datasets using bio-omics tools and statistical analysis.



### Motivations

The motivation behind this project was skepticism about the PandaOmics weighting algorithm and how it combines datasets when default parameters are set in the application. Because of these non-transparent practices, concerns were raised about potentially conflicting datasets for our research. This independently built project was used to verify findings using the raw, identified datasets and statistical methods that could be reproduced.



### Pipeline

The MS analysis Python pipeline (ms\_analysis.py) handles the workflow as follows;

* Loading raw expression, metadata CSV's downloaded from each GEO dataset
* Determine if a dataset is pre or post-log transformed, to avoid issues from double-log-transforming a dataset
* Map Ensembl gene IDs to HGNC using Biomart Ensembl Genes 115
* Separate samples into MS or health control using pandas per dataset
* Compute the log2FC transformation conditionally for each dataset's transformation dataset
* Run an independent t-test per gene per data set (individually) to avoid any batch effects
* Summarize upregulation/downregulation counts per the gene and tissue type into a CSV
* Visualize plots of log2FC in datasets and up/down regulation



### Libraries used

pandas, numpy, scipy, seaborn, matplotlib



## Results

* 90% of the 10 raw GEO datasets successfully came out of the pipeline; One was unsuccessful due to, upon further analysis of the raw file, being identified as a murine model
* Concerning Calprotectin (S100A8/9), 7 out of 9 datasets, or 77%, confirmed significant upregulation in MS patients versus controls. A binomial meta-analysis concluded p to be statistically significant at 0.09. While not conventionally significant, these small-scale prompts further investigation, particularly due to the magnitude of the Log2 Fold Change across both CNS lesions and peripheral blood samples.
* Concerning ALOX15B, 6 out of the 9 datasets, or 66%, confirmed upregulation.



|File|Description|
|-|-|
|results.csv|Each dataset's results including the log2FC, means, sample sizes and p-values.|
|summary.csv|Summary, aggregated with gene and tissue type.|
|plot1\_log2fc\_per\_dataset.png|Log2fold change per dataset per gene, colored by tissue type|
|plot2\_updown\_summary.png|Upregulation vs downregulation counts by tissue|

|poster.pdf|The presentable poster that was made for this research.|
|-|-|





## References

* MS Global Impact: Multiple Sclerosis International Federation. *Atlas of MS 3rd Edition.* 2020
* Reich, Daniel S., et al. "Multiple Sclerosis." *NEJM*, vol. 378, no. 2, 2018, pp. 169–80
* Safiadeh, Elham, et al. "Arachidonate 15-Lipoxygenase Type B Is Highly Expressed in Active MS Lesions." *Journal of Neuroinflammation*, vol. 15, no. 1, 2018
* Foell, Dirk, et al. "S100 Proteins in Health and Disease." *Journal of Leukocyte Biology*, vol. 83, no. 4, 2008
* Insilico Medicine. PandaOmics Platform. 2026
* NCBI Gene Expression Omnibus (GEO). 2026

