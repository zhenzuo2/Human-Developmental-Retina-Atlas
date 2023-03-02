MODEL_H5 = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/earlyPRPC/models"
REGIONS = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/earlyPRPC/auxiliary/filtered.nonpeaks.bed"
GENOME = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
CHR_FOLD_PATH = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json"
OUTPUT_PREFIX = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_/chrombpnet_marginal_footprinting/earlyPRPC/"
MOTIFS_TO_PWM = ""
chrombpnet footprints -m MODEL_H5 -r REGIONS -g GENOME -fl CHR_FOLD_PATH -op OUTPUT_PREFIX -pwm_f MOTIFS_TO_PWM
                             [-bs BATCH_SIZE] [--ylim YLIM]