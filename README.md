# scNEAT
scNEAT is a pipeline to analyze NEAT-seq data containing five analyses:  
1) Clustering NEAT-seq data by using scHCAE  
2) Target prediction of TFs by using scMF 
3) Cell-type specific DE analysis  
4) Regression analysis to detect confident target genes of TFs
5) Regression analysis to detect co-target genes of two TFs.

# Dependencies
Python:  
Python 3.8.1  
Pytorch 1.6.0  
Scanpy 1.6.0  
SKlearn 0.22.1  
Numpy 1.18.1  
h5py 2.9.0  

#R:  
Seurat 4.2.0  
rhdf5 2.38.1   
chromVAR 1.16.0  
JASPAR2016 1.22.0  
chromVARmotifs 0.2.0  
motifmatchr 1.16.0
SingleCellExperiment 1.16.0
BiocParallel 1.28.3
