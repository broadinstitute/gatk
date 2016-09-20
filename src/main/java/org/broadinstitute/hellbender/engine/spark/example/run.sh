#!/bin/bash

if [ -z "$SPARK_HOME" ]; then
  echo "SPARK_HOME not set"
  exit
fi

CLASS=com.intel.genomicsdb.GenomicsDBSparkFactory
GENOMICSDB_SPARK_JAR=genomicsdb-sparkapi-0.1.jar
GENOMICSDB_SPARK_HOME="$(cd `dirname $0`; pwd)"
GENOMICSDB_SPARK_JAR=genomicsdb-sparkapi-0.1.jar
GENOMICSDB_BASE_JAR=$GENOMICSDB_SPARK_HOME/../../bin/genomicsdb.jar
HTSJDK_JAR=$HOME/Downloads/htsjdk-2.4.1.jar
BREEZE_JAR=$HOME/Downloads/breeze-math_2.10-0.4.jar:$HOME/Downloads/breeze-core_2.10-0.2.jar

CLASSPATH=$GENOMICSDB_SPARK_HOME/../target/$GENOMICSDB_SPARK_JAR
CLASSPATH+=:$HTSJDK_JAR:$GENOMICSDB_BASE_JAR:$BREEZE_JAR

$SPARK_HOME/bin/spark-submit \
  --deploy-mode client \
  --class $CLASS \
  --driver-class-path $CLASSPATH \
  --driver-library-path $LD_LIBRARY_PATH \
  --conf spark.executor.extraLibraryPath=$LD_LIBRARY_PATH \
  --conf spark.executor.extraClassPath=$CLASSPATH \
  --conf spark.tiledb.workspace.dir=$TILEDB_WORKSPACE \
  target/$GENOMICSDB_SPARK_JAR \
  "$@"
