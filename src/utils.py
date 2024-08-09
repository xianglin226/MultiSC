from munkres import Munkres
import os
import sys
import math
import scanpy as sc
from scipy import stats, spatial, sparse
from scipy.linalg import norm
import scipy
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import random
import torch
import torch.nn as nn
import torch.nn.init as init
import torch.utils.data as data
from sklearn.neighbors import kneighbors_graph

# def cluster_acc(y_true, y_pred):
    # """
    # Calculate clustering accuracy. Require scikit-learn installed
    # # Arguments
        # y: true labels, numpy.array with shape `(n_samples,)`
        # y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    # # Return
        # accuracy, in [0,1]
    # """
    # y_true = y_true.astype(np.int64)
    # assert y_pred.size == y_true.size
    # D = max(y_pred.max(), y_true.max()) + 1
    # w = np.zeros((D, D), dtype=np.int64)
    # for i in range(y_pred.size):
        # w[y_pred[i], y_true[i]] += 1
    # from sklearn.utils.linear_assignment_ import linear_assignment
    # ind = linear_assignment(w.max() - w)
    # return sum([w[i, j] for i, j in ind]) * 1.0 / y_pred.size

def GetCluster(X, res, n):
    adata0=sc.AnnData(X)
    if adata0.shape[0]>200000:
       np.random.seed(adata0.shape[0])#set seed 
       adata0=adata0[np.random.choice(adata0.shape[0],200000,replace=False)] 
    sc.pp.neighbors(adata0, n_neighbors=n, use_rep="X")
    sc.tl.louvain(adata0,resolution=res)
    Y_pred_init=adata0.obs['louvain']
    Y_pred_init=np.asarray(Y_pred_init,dtype=int)
    if np.unique(Y_pred_init).shape[0]<=1:
        #avoid only a cluster
        exit("Error: There is only a cluster detected. The resolution:"+str(res)+"is too small, please choose a larger resolution!!")
    else: 
        print("Estimated n_clusters is: ", np.shape(np.unique(Y_pred_init))[0]) 
    return(np.shape(np.unique(Y_pred_init))[0])

def torch_PCA(X, k, center=True, scale=False):
    X = X.t()
    n,p = X.size()
    ones = torch.ones(n).cuda().view([n,1])
    h = ((1/n) * torch.mm(ones, ones.t())) if center else torch.zeros(n*n).view([n,n])
    H = torch.eye(n).cuda() - h
    X_center =  torch.mm(H.double(), X.double())
    covariance = 1/(n-1) * torch.mm(X_center.t(), X_center).view(p,p)
    scaling = torch.sqrt(1/torch.diag(covariance)).double() if scale else torch.ones(p).cuda().double()
    scaled_covariance = torch.mm(torch.diag(scaling).view(p,p), covariance)
    eigenvalues, eigenvectors = torch.eig(scaled_covariance, True)
    components = (eigenvectors[:, :k])
    #explained_variance = eigenvalues[:k, 0]
    return components
    
def best_map(L1,L2):
            #L1 should be the groundtruth labels and L2 should be the clustering labels we got
            Label1 = np.unique(L1)
            nClass1 = len(Label1)
            Label2 = np.unique(L2)
            nClass2 = len(Label2)
            nClass = np.maximum(nClass1,nClass2)
            G = np.zeros((nClass,nClass))
            for i in range(nClass1):
                ind_cla1 = L1 == Label1[i]
                ind_cla1 = ind_cla1.astype(float)
                for j in range(nClass2):
                    ind_cla2 = L2 == Label2[j]
                    ind_cla2 = ind_cla2.astype(float)
                    G[i,j] = np.sum(ind_cla2 * ind_cla1)
            m = Munkres()
            index = m.compute(-G.T)
            index = np.array(index)
            c = index[:,1]
            newL2 = np.zeros(L2.shape)
            for i in range(nClass2):
                newL2[L2 == Label2[i]] = Label1[c[i]]
            return newL2
 
def geneSelection(data, threshold=0, atleast=10, 
                  yoffset=.02, xoffset=5, decay=1.5, n=None, 
                  plot=True, markers=None, genes=None, figsize=(6,3.5),
                  markeroffsets=None, labelsize=10, alpha=1, verbose=1):
    
    if sparse.issparse(data):
        zeroRate = 1 - np.squeeze(np.array((data>threshold).mean(axis=0)))
        A = data.multiply(data>threshold)
        A.data = np.log2(A.data)
        meanExpr = np.zeros_like(zeroRate) * np.nan
        detected = zeroRate < 1
        meanExpr[detected] = np.squeeze(np.array(A[:,detected].mean(axis=0))) / (1-zeroRate[detected])
    else:
        zeroRate = 1 - np.mean(data>threshold, axis=0)
        meanExpr = np.zeros_like(zeroRate) * np.nan
        detected = zeroRate < 1
        mask = data[:,detected]>threshold
        logs = np.zeros_like(data[:,detected]) * np.nan
        logs[mask] = np.log2(data[:,detected][mask])
        meanExpr[detected] = np.nanmean(logs, axis=0)

    lowDetection = np.array(np.sum(data>threshold, axis=0)).squeeze() < atleast
    zeroRate[lowDetection] = np.nan
    meanExpr[lowDetection] = np.nan
            
    if n is not None:
        up = 10
        low = 0
        for t in range(100):
            nonan = ~np.isnan(zeroRate)
            selected = np.zeros_like(zeroRate).astype(bool)
            selected[nonan] = zeroRate[nonan] > np.exp(-decay*(meanExpr[nonan] - xoffset)) + yoffset
            if np.sum(selected) == n:
                break
            elif np.sum(selected) < n:
                up = xoffset
                xoffset = (xoffset + low)/2
            else:
                low = xoffset
                xoffset = (xoffset + up)/2
        if verbose>0:
            print('Chosen offset: {:.2f}'.format(xoffset))
    else:
        nonan = ~np.isnan(zeroRate)
        selected = np.zeros_like(zeroRate).astype(bool)
        selected[nonan] = zeroRate[nonan] > np.exp(-decay*(meanExpr[nonan] - xoffset)) + yoffset
                
    if plot:
        if figsize is not None:
            plt.figure(figsize=figsize)
        plt.ylim([0, 1])
        if threshold>0:
            plt.xlim([np.log2(threshold), np.ceil(np.nanmax(meanExpr))])
        else:
            plt.xlim([0, np.ceil(np.nanmax(meanExpr))])
        x = np.arange(plt.xlim()[0], plt.xlim()[1]+.1,.1)
        y = np.exp(-decay*(x - xoffset)) + yoffset
        if decay==1:
            plt.text(.4, 0.2, '{} genes selected\ny = exp(-x+{:.2f})+{:.2f}'.format(np.sum(selected),xoffset, yoffset), 
                     color='k', fontsize=labelsize, transform=plt.gca().transAxes)
        else:
            plt.text(.4, 0.2, '{} genes selected\ny = exp(-{:.1f}*(x-{:.2f}))+{:.2f}'.format(np.sum(selected),decay,xoffset, yoffset), 
                     color='k', fontsize=labelsize, transform=plt.gca().transAxes)

        plt.plot(x, y, color=sns.color_palette()[1], linewidth=2)
        xy = np.concatenate((np.concatenate((x[:,None],y[:,None]),axis=1), np.array([[plt.xlim()[1], 1]])))
        t = plt.matplotlib.patches.Polygon(xy, color=sns.color_palette()[1], alpha=.4)
        plt.gca().add_patch(t)
        
        plt.scatter(meanExpr, zeroRate, s=1, alpha=alpha, rasterized=True)
        if threshold==0:
            plt.xlabel('Mean log2 nonzero expression')
            plt.ylabel('Frequency of zero expression')
        else:
            plt.xlabel('Mean log2 nonzero expression')
            plt.ylabel('Frequency of near-zero expression')
        plt.tight_layout()
        
        if markers is not None and genes is not None:
            if markeroffsets is None:
                markeroffsets = [(0, 0) for g in markers]
            for num,g in enumerate(markers):
                i = np.where(genes==g)[0]
                plt.scatter(meanExpr[i], zeroRate[i], s=10, color='k')
                dx, dy = markeroffsets[num]
                plt.text(meanExpr[i]+dx+.1, zeroRate[i]+dy, g, color='k', fontsize=labelsize)
    
    return selected


def generate_random_pair(adt, label_cell_indx, low, high, num1, num2):
    """
    Generate random pairwise constraints.
    """
    ml_ind1, ml_ind2 = [], []
    cl_ind1, cl_ind2 = [], []
    adt = np.array(adt)

    def check_ind(ind1, ind2, ind_list1, ind_list2):
        for (l1, l2) in zip(ind_list1, ind_list2):
                if ind1 == l1 and ind2 == l2:
                    return True
        return False

    lim = 1000000
    while (num1 > 0 or num2 >0) and lim >0:
        tmp1 = random.choice(label_cell_indx)
        tmp2 = random.choice(label_cell_indx)
        if tmp1 == tmp2:
            continue
        if check_ind(tmp1, tmp2, ml_ind1, ml_ind2):
            continue
        if adt[tmp1] > high and adt[tmp2] > high:
          if num1 >0:
             ml_ind1.append(tmp1)
             ml_ind2.append(tmp2)
             num1 -= 1
        elif adt[tmp1] < low and adt[tmp2] < low:
          if num1 >0:
             ml_ind1.append(tmp1)
             ml_ind2.append(tmp2)
             num1 -= 1
        elif adt[tmp1] > high and adt[tmp2] < low:
          if num2 >0:
             cl_ind1.append(tmp1)
             cl_ind2.append(tmp2)
             num2 -= 1
        elif adt[tmp2] > high and adt[tmp1] < low:
          if num2 >0:
             cl_ind1.append(tmp1)
             cl_ind2.append(tmp2)
             num2 -= 1
        else:
          continue
          lim-=1
               
    ml_ind1, ml_ind2 = np.array(ml_ind1), np.array(ml_ind2)
    cl_ind1, cl_ind2 = np.array(cl_ind1), np.array(cl_ind2)
    ml_index = np.random.permutation(ml_ind1.shape[0])
    cl_index = np.random.permutation(cl_ind1.shape[0])
    ml_ind1 = ml_ind1[ml_index]
    ml_ind2 = ml_ind2[ml_index]
    cl_ind1 = cl_ind1[cl_index]
    cl_ind2 = cl_ind2[cl_index]
    print("Number of ML: ",ml_ind1.shape)
    print("Number of CL: ",cl_ind1.shape)
    return ml_ind1, ml_ind2, cl_ind1, cl_ind2

def clr_normalize_each_cell(adata):
    """Normalize count vector for each cell, i.e. for each row of .X"""

    def seurat_clr(x):
        # TODO: support sparseness
        s = np.sum(np.log1p(x[x > 0]))
        exp = np.exp(s / len(x))
        return np.log1p(x / exp)
    
    adata.raw = adata.copy()
    
    adata.X = adata.X.astype(float)
    sc.pp.normalize_per_cell(adata)
    adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)

    # apply to dense or sparse matrix, along axis. returns dense matrix
    adata.X = np.apply_along_axis(
        seurat_clr, 1, (adata.raw.X.A if scipy.sparse.issparse(adata.raw.X) else adata.raw.X)
    )
    return adata.X
