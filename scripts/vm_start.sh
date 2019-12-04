#!/usr/bin/env bash
VM=${1:-sam-gpu2}
shift 1
ZONE=us-central1-a
MAX_TRIES=1000
COUNTER=0
while [[ $COUNTER -lt $(( $MAX_TRIES )) ]]; do
    gcloud compute instances start $VM --zone $ZONE
    if [[ $? -eq 0 ]]
    then
      echo "Successfully started vm: ${VM} after ${COUNTER} attempts."
      break
    else
      let COUNTER=COUNTER+1
      echo "Could not start vm: ${VM}, unsuccessful attempt: ${COUNTER}."
    fi
done

gcloud compute ssh $VM --zone $ZONE


