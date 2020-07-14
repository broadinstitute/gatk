WORKING_DIR=/home/travis/build/broadinstitute

docker pull ga4gh/htsget-ref:1.1.0
docker container run -d --name htsget-server -p 3000:3000 --env HTSGET_PORT=3000 --env HTSGET_HOST=http://0.0.0.0:3000 \
      -v $WORKING_DIR/gatk/src/test/resources/org/broadinstitute/hellbender/engine:/data \
      ga4gh/htsget-ref:1.1.0 \
      ./htsref -config /data/htsget_config.json
docker container ls -a

curl http://localhost:3000