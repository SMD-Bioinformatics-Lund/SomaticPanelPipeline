DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PIPELINE_DEST="/fs1/viktor/SomaticPanelPipeline_dsl2/"
DEST_HOST="rs-fs1.lunarc.lu.se"


# Copy pipeline script
scp -r $DIR/workflows $DEST_HOST:$PIPELINE_DEST
scp -r $DIR/subworkflows $DEST_HOST:$PIPELINE_DEST
scp -r $DIR/modules $DEST_HOST:$PIPELINE_DEST


# Copy configuration file
scp $DIR/configs/nextflow.hopper.config $DEST_HOST:$PIPELINE_DEST/nextflow.config

# Copy other files
scp -r $DIR/bin $DEST_HOST:$PIPELINE_DEST

scp $DIR/scenario.yaml $DEST_HOST:/fs1/resources/ref/hg38/annotation_dbs/varlociraptor/scenario_ffpe.yaml
git rev-parse HEAD > git.hash
scp $DIR/git.hash $DEST_HOST:$PIPELINE_DEST
