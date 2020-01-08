#!/bin/bash

export REPO_PATH=/staging/as/kaliappa/tmp/rb_utils

snakemake \
--snakefile $REPO_PATH/snakefiles/pe_cnv/Snakefile \
--configfile $REPO_PATH/snakefiles/pe_cnv/config.yaml \
--printshellcmds \
--keep-going \
--rerun-incomplete \
--cluster-config $REPO_PATH/snakefiles/pe_cnv/cluster.yaml \
--cores 10 \
--cluster 'sbatch --partition={cluster.partition} --ntasks={cluster.cores} --mem={cluster.mem} --time={cluster.time} -o {cluster.logout} -e {cluster.logerror}' \
$@
