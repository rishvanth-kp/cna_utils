#!/bin/bash

export REPO_PATH=/scratch2/kaliappa/cnv/cna_utils/snakefiles/cnv

if [ ! -e "slurm_logs" ]
then
  mkdir "slurm_logs"
fi

snakemake \
--snakefile $REPO_PATH/Snakefile \
--configfile $REPO_PATH/config.yaml \
--printshellcmds \
--keep-going \
--rerun-incomplete \
--cluster-config $REPO_PATH/cluster.yaml \
--cores 50 \
--cluster 'sbatch --ntasks={cluster.tasks} --cpus-per-task={cluster.cores} --mem={cluster.mem} --time={cluster.time} --account={cluster.acc} --output {cluster.logout} --error {cluster.logerror}' \
$@
