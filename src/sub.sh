#!/bin/bash -l
#SBATCH -J 1
#SBATCH -p gpu 
#SBATCH --gres=gpu:1
#SBATCH --qos=low                         
#SBATCH --time 72:00:00

module load CUDA/11.4.1
 conda activate scbert
f=/project/zhiwei/sj225/MultiSC/Input/GSE178707_neatseq_lane2.h5

python -u run_scMultiCluster3.py\
 --n_clusters 6\
 --data_file $f\
 --save_dir /project/zhiwei/sj225/MultiSC/result3/line2/FOXP3/2000_3000_ml1/\




