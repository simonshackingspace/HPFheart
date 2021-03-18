#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate schpf_p37

DATA_DIR=/home/zz2565/data/thesis_data/diff_seqs/
SCRNA_FILE=${DATA_DIR}sc_data_filtered.csv
PRETRAIN_DIR=${DATA_DIR}scHPF_pretrain/
TRAIN_DIR=${DATA_DIR}scHPF_train/
SCORE_DIR=${DATA_DIR}scHPF_score/
BEST_SCORE_DIR=thesis_output/scHPF_best_score/
bestK=0

# umi file for scHPF
python mkgl.py $SCRNA_FILE

# scHPF preprocessing
scHPF prep -i ${SCRNA_FILE}.umi -o $PRETRAIN_DIR -m 10

while true; do
	echo "input number of factors: "	
	read K 
	rm -r $TRAIN_DIR
	rm -r $SCORE_DIR
	if [ $K -eq 0 ]; then
		echo "best K=${bestK} chosen"
		break
	fi
	scHPF train -i ${PRETRAIN_DIR}filtered.mtx -o $TRAIN_DIR -k $K -e 0.01
	scHPF score -m ${TRAIN_DIR}*joblib -o $SCORE_DIR -g ${PRETRAIN_DIR}genes.txt
	Rscript schpf_finding_K.R ${SCORE_DIR}cell_score.txt $K ${DATA_DIR}celltypes.csv
	bestK=$K
done
scHPF train -i ${PRETRAIN_DIR}filtered.mtx -o $TRAIN_DIR -k $bestK -t 5
scHPF score -m ${TRAIN_DIR}*joblib -o ${BEST_SCORE_DIR} -g ${PRETRAIN_DIR}genes.txt
rm -r ${TRAIN_DIR}
Rscript schpf_finding_K.R ${BEST_SCORE_DIR}cell_score.txt $bestK ${DATA_DIR}celltypes.csv
	
conda deactivate
