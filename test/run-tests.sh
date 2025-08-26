 #!/usr/bin/env bash

set -xeuo pipefail

script_dir="$(dirname $0)"
DOCKER=${1:-no}

if [ "$DOCKER" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-profile gh'
else
    export SINGULARITY_FAKEROOT=1
    docker_flag='-profile local -with-singularity'
fi

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    -c "$script_dir"/nextflow.config \
    --sample_sheet "$script_dir"/inputs/sample-sheet.csv \
    --inputs "$script_dir"/inputs \
    --outputs "$script_dir"/outputs
