# Docker Instructions

## Build docker image
```shell
docker build . -t leafcutter2
```

## Check docker image
```shell
docker run leafcutter2 --help
```

## Run docker container
```shell
docker run \
    -v ./example:/workspace/junctions \
    -v ./example/data:/workspace/data \
    -v ./example/annotation:/workspace/annotation \
    -v ./output:/workspace/output \
    leafcutter2 \
    -j junctions/junction_files.txt \
    -r output \
    -A annotation/chr10.gtf.gz \
    -G annotation/chr10.fa.gz
```
