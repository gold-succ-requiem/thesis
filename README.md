# In silico network-based drug repositioning through similarity fusion and adverse drug event penalisation

## Background
A network-based drug repositioning system that incorporates drug structure, drug-induced pathway, target sequence, and target function similarity data. A key aspect of the project is to incorporate side effect and toxicity data. Rather than recommending repositioning candidates based on the similarity of their side effects, it penalises candidates with more severe side effects. Side effect frequency is used as a proxy for severity. It is assumed side effects with higher frequencies are less severe than those with lower frequencies.

## Features

### At a glance
- similarity aggregation via similarity network fusion
- clustering analysis via spectral clustering
- reprioritisation of recommended drugs based on side effects

### Aggregation by similarity network fusion
#### Input
Similarity matrices for drug and disease data, plus corresponding drug IDs in the RepoDB dataset.

- `chem_similarity.csv`: drug structure similarity
- `target_similarity.csv`: drug target tertiary sequence similarity
- `path_similarity_v2.csv`: drug-induced pathway similarity
- `GO_Sim_BP_combined.csv`: drug target GO term (biological process) similarity
- `repoDB_full.csv`: corresponding drug IDs in RepoDB dataset

#### Processing
- `snf2.r`

Performs the following functions:

- **Similarity matrix preprocessing.** Sets matrix rownames; resizes matrix dimensions to intersection of drug IDs in matrix and RepoDB; filters rows and columns containing entirely NA values.
- **Combination.** Generates list of combinations of matrices (e.g. structure only; structure and target; structure, pathway, and BP; all four). Gives

![2^n  - 1](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

combinations of matrices. In our case, n = 4.

- **Similarity network fusion.** Based on [Wang et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/24464287/). d
- **Fusion matrix filtering.** 

#### Output
- `W.combn.x.rds`
- `W.x.rds`

Refer to the following table for x values:

x|Quantile
-|-
1|0th
3|50th
4|75th

Thus, `x` maps to quantile at which fused matrices were filtered. e.g. `W.3.rds` will contain a list of fusion matrices where similarities less than the median have been converted to 0.

### Reprioritisation based on side effects

# References
