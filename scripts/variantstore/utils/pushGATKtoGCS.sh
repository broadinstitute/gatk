JAR=$(ls -t build/libs/*-SNAPSHOT-local.jar | head -1)
DT=$(date '+%Y%m%d')

branch=$(git symbolic-ref HEAD 2>/dev/null)
branch=${branch#refs/heads/} 

DEST="gs://gvs-internal/bigquery-jointcalling/jars/${branch}_${DT}/"

gsutil cp $JAR $DEST

echo "Copied to $DEST$(basename $JAR)"
