## BINRES: Bayesian integrative region segmentation in spatially resolved transcriptomic studies

The R package BINRES is built to implement the nonparametric Bayesian method named BINRES to carry out the region segmentation for a tissue section by integrating all the three types of data generated during the study --- gene expressions, spatial coordinates, and the histology image. BINRES captures the spatial dependence of neighboring spots and does not require a prespecified region number. It also combines the image and the gene expressions whose contribution weights can be flexibly adjusted in a data-adaptive manner. The computationally scalable extension BINRES-fast is developed for large-scale studies. In this package, a partially collapsed Gibbs sampler is carefully designed for Bayesian posterior inference. BINRES can be installed in Windows, Linux, and Mac OS.

For technical details, please refer to our paper currently accepted in *Journal of the American Statistical Association*: Yinqiao Yan and Xiangyu Luo, "Bayesian integrative region segmentation in spatially resolved transcriptomic studies" with DOI: xxx and URL: [xxx](xxx). 

The code for reproducibility in the paper can be downloaded through [xxx](xxx).



## Prerequisites and Installation

1. R version >= 4.1.3.
2. CRAN package: mvtnorm (>=1.1.3), MCMCpack (>=1.6.3), rBeta2009 (>=1.0), truncnorm (>=1.0.8), stats (>=4.1.3)
3. Install the package BINRES.

```R
devtools::install_github("yinqiaoyan/BINRES")
```

## Example Code

The following code shows an example (Simulation I in the paper) that runs the main function "BINRES" in our package.

```R
library(BINRES)
library(aricode)
library(ggplot2)
# Import example data
# (1) coord: Spatial coordinates
# (2) image_data: Simulated image data (normally distributed)
# (3) gene_data: Simulated gene expression data (excess zeros and log-normally distributed)
# (4) true_label: True region labels of all spots
# (5) marker_id: ID number of marker genes
data(example_data_BINRES)
# Dimension of spatial coordinates
dim(coord)
# Dimension of image data
dim(image_data)
# Dimension of gene expression data
dim(gene_data)
# Auxiliary functions
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}
# --- run BINRES ---
# Total execution time is about 4.5 minutes 
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
res_list = BINRES(image_data = image_data, gene_data = gene_data, coord = coord,
                  platform="ST", num_init=5,
                  Is_logNormPcaGene=FALSE,
                  numOfMCMC=1000, burnIn=500, print_gap=10,
                  Is_warm_start=FALSE, Is_random_seed=TRUE, random_seed=30)
# Execution time
res_list$exeTime
# Posterior mode of consensus clustering C and marker-gene indicator \gamma
clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)
gamma_mode = apply(res_list$gamma_g_mcmc, 2, getmode)
# 95% credible interval for spatial interaction parameter \beta
quantile(res_list$pottsBeta_mcmc, c(0.025, 0.975))
# True marker gene ID
cat("True marker gene ID:", marker_id)
# Estimated marker gene ID
cat("Estimated marker gene ID:", which(gamma_mode == 1))
# Compared with true labels
table(clIds_mode, true_label)
cat("ARI value:", ARI(clIds_mode, true_label))
# Visualization
tmpc = clIds_mode
tmpc2 = tmpc
tmpc2[tmpc == 1] = 2
tmpc2[tmpc == 2] = 4
tmpc2[tmpc == 3] = 3
tmpc2[tmpc == 4] = 1
tmpc = tmpc2
plot_color=c("#CC6677", "#53B8EA", "#E2C845", "#03AF3C")
plot_dat <- data.frame(x = coord[,1], y = coord[,2], c = tmpc)
p <- ggplot(data = plot_dat, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 3) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values=plot_color)
par(ask=FALSE)
print(p)
```

The following code shows another example (Mouse coronal brain section data in the paper) that runs the main function "BINRES_fast" in our package.

```R
library(BINRES)
library(aricode)
library(ggplot2)
# Import example data
# (1) coord: Spatial coordinates
# (2) image_data: preprocessed image data
# (3) gene_data_pc: preprocessed gene expression matrix
# obtained by normalizing ST raw count matrix, taking logarithm, and conducting PCA
# (4) true_label: Biotype annotations (true region labels) of all spots
data(example_data_BINRES_fast)
# Dimension of spatial coordinates
dim(coord)
# Dimension of image data
dim(image_data)
# Dimension of gene expression data
dim(gene_data_pc)
# Auxiliary functions
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}
# --- run BINRES ---
# Total execution time is about 30 minutes 
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
res_list_fast = BINRES_fast(image_data = image_data, gene_data_pc = gene_data_pc, coord = coord,
                            platform="Visium", num_init=15, 
                            numOfMCMC=3000, burnIn=1500, print_gap=50,
                            Is_warm_start=TRUE, Is_random_seed=TRUE, random_seed=78)
# Execution time
res_list_fast$exeTime
# Posterior mode of consensus clustering C
clIds_mode = apply(res_list_fast$clIds_mcmc, 2, getmode)
# 95% credible interval for spatial interaction parameter \beta
quantile(res_list_fast$pottsBeta_mcmc, c(0.025, 0.975))
# Compared with true labels
table(clIds_mode, true_label)
cat("ARI value:", ARI(clIds_mode, true_label))
# Visualization

```

or you can simply run

```R
library(BINRES)
example("BINRES")
example("BINRES_fast")
```

## Remarks

- If you have any questions regarding this package, please contact Yinqiao Yan at [yanyinqiao@ruc.edu.cn](mailto:yanyinqiao@ruc.edu.cn).