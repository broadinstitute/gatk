# this script starts a docker container running a GA4GH reference htsget server for use in htsget tests
WORKING_DIR=/home/travis/build/broadinstitute

docker pull ga4gh/htsget-ref:1.1.0
docker container run -d --name htsget-server -p 3000:3000 --env HTSGET_PORT=3000 --env HTSGET_HOST=http://127.0.0.1:3000 \
      -v $WORKING_DIR/gatk/src/test/resources/org/broadinstitute/hellbender/engine:/data \
      ga4gh/htsget-ref:1.1.0 \
      ./htsref -config /data/htsget_config.json
docker container ls -a

curl http://localhost:3000
