#!/usr/bin/env bash
VM=${1:-sam-gpu2}
shift 1
ZONE=us-central1-a
MAX_TRIES=1000
COUNTER=0
while [[ $COUNTER -lt $(( $MAX_TRIES )) ]]; do
    sleep 1s
    gcloud compute instances start $VM --zone $ZONE
    if [[ $? -eq 0 ]]
    then
      echo "Potentially started vm: ${VM} after ${COUNTER} attempts."
      gcloud compute ssh $VM --zone $ZONE
      if [[ $? -eq 0 ]]
      then
        break
      else
        let COUNTER=COUNTER+1
        echo "Actually, no. Could not start vm: ${VM}, unsuccessful attempt: ${COUNTER}."
        sleep 1s
      fi
    else
      let COUNTER=COUNTER+1
      sleep 1s
      echo "Could not start vm: ${VM}, unsuccessful attempt: ${COUNTER}."
      sleep 1s
    fi
done




