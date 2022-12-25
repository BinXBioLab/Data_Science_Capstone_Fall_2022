PIPELINE_SCRIPT=/pipeline/python/pipeline.py
INTER=/pipeline/inter

echo "Running Preprocessing"
python $ROOT/python/preprocessing.py \
    --in_dir=/pipeline/raw_data \
    --inter_dir=$INTER

# Work in progress
echo "Running Main Pipeline"
python $PIPELINE_SCRIPT run-liger \
    --rscript=/pipeline/r/run_liger.R \
    --counts=$INTER/pasca_log1p.parquet \
    --metadata=$INTER/pasca_log1p.metadata.csv \
    --output=$INTER/pasca.liger.parquet

python $PIPELINE_SCRIPT post-liger-pre-singler \
    --anndata=$INTER/pasca_log1p.h5ad \
    --liger=$INTER/pasca.liger.parquet \
    --output=$INTER/pasca.cpm.parquet

python $PIPELINE_SCRIPT run-singler \
    --rscript=/pipeline/r/run_singler.R \
    --output=$INTER/pasca_nowakowski_med_2022-09-25_singler.csv

python $PIPELINE_SCRIPT post-singler-pre-mast \
    --anndata=$INTER/pasca.cpm.parquet \
    --annotation=$INTER/pasca_nowakowski_med_2022-09-25_singler.csv \
    --inter=$INTER

python $PIPELINE_SCRIPT run-mast \
    --rscript=/pipeline/r/run_mast.R \
    --anndata=$INTER/???.h5ad \
    --metadata=$INTER/nowakowski_EN.mast_data.csv \
    --output=$INTER/nowakowski_EN.mast_ent_de.csv

python $PIPELINE_SCRIPT post-mast \
    --ent_de=$INTER/nowakowski_EN.mast_ent_de.csv

