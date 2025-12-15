#!/bin/bash
#SBATCH --gres=gpu:3
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -J L_trnew
#SBATCH -p gpu
#SBATCH -o /mnt/scratch/jugipalace/foldseek_new_3di/log_3di_new/L_trnew.%j.out
#SBATCH -e /mnt/scratch/jugipalace/foldseek_new_3di/log_3di_new/L_trnew.%j.err

cd /mnt/scratch/jugipalace/foldseek_new_3di/training_3di_gpu_new
source ~/miniconda3/etc/profile.d/conda.sh
conda activate fold3di
./learnAlphabet.sh 20 100 data/pdbs_train.txt data/pdbs_val.txt tmp/
