time_gatk() {
  GATK_ARGS=$1
  NUM_EXECUTORS=$2
  EXECUTOR_CORES=$3
  EXECUTOR_MEMORY=$4
  DRIVER_MEMORY=$5
  if [ -n "$GCS_CLUSTER" ]; then
    API_KEY_ARGS="-apiKey $API_KEY"
    SPARK_RUNNER_ARGS="--sparkRunner GCS --cluster $GCS_CLUSTER"
  else
    API_KEY_ARGS=""
    SPARK_RUNNER_ARGS="--sparkRunner SPARK --sparkMaster yarn-client --sparkSubmitCommand spark2-submit"
  fi
  COMMAND=$(echo $GATK_ARGS | awk '{print $1}')
  RESULTS_CSV=results/$(basename "$0" .sh).csv
  mkdir -p results
  LOG=logs/${COMMAND}_$(date +%Y%m%d_%H%M%S).log
  mkdir -p logs
  echo "$GATK_HOME/gatk-launch $GATK_ARGS $API_KEY_ARGS \\
    -- \\
    $SPARK_RUNNER_ARGS \\
    --num-executors $NUM_EXECUTORS --executor-cores $EXECUTOR_CORES --executor-memory $EXECUTOR_MEMORY \\
    --driver-memory $DRIVER_MEMORY \\
    --conf spark.dynamicAllocation.enabled=false" \
  > $LOG
  $GATK_HOME/gatk-launch $GATK_ARGS $API_KEY_ARGS \
    -- \
    $SPARK_RUNNER_ARGS \
    --num-executors $NUM_EXECUTORS --executor-cores $EXECUTOR_CORES --executor-memory $EXECUTOR_MEMORY \
    --driver-memory $DRIVER_MEMORY \
    --conf spark.dynamicAllocation.enabled=false \
  >> $LOG 2>&1
  RC=$?
  DURATION_MINS=$(grep 'Elapsed time' $LOG | grep -Eow "[0-9]+\.[0-9][0-9]")
  if [ ! -e $RESULTS_CSV ]; then
    echo 'Command,Number of executors,Executor cores,Executor memory,Exit code,Time (mins)' > $RESULTS_CSV
  fi
  echo "$GATK_ARGS,$NUM_EXECUTORS,$EXECUTOR_CORES,$EXECUTOR_MEMORY,$RC,$DURATION_MINS" >> $RESULTS_CSV
}
