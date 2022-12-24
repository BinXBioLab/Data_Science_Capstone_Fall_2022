ROOT=/pipeline

# echo "Running Preprocessing"
# python $ROOT/python/preprocessing.py --in_dir=$ROOT/raw_data --inter_dir=$ROOT/inter

# Work in progress
echo "Running Main Pipeline"
python $ROOT/python/pipeline.py run-liger
python $ROOT/python/pipeline.py post-liger-pre-singler
python $ROOT/python/pipeline.py run-singler
python $ROOT/python/pipeline.py post-singler-pre-mast
python $ROOT/python/pipeline.py run-mast
python $ROOT/python/pipeline.py post-mast

