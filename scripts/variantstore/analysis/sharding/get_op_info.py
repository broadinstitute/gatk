from dateutil.parser import parse
import sys
from pprint import pprint

from googleapiclient.discovery import build
from oauth2client.client import GoogleCredentials

credentials = GoogleCredentials.get_application_default()
service = build('lifesciences', 'v2beta', credentials=credentials)

for opid in sys.stdin:
    opid = opid.strip();

    request = service.projects().locations().operations().get(name=opid)
    response = request.execute()

    start_ts = None
    end_ts = None

    m = response['metadata']
    for e in m['events']:
        timestamp = parse(e['timestamp'])
        desc = e['description']
        
        # look for latest "assigned" and latest "Worker released" as the timestamps.. have to check timestamps...
        if ("assigned" in desc and (start_ts is None or timestamp >= start_ts)):
            start_ts = timestamp
            
        if ("Worker released" in desc and (end_ts is None or timestamp >= end_ts)):
            end_ts = timestamp

    if (start_ts is not None):
        print(f"{opid}\t{start_ts}\t{end_ts}")
    else:
        print(f"{opid}")
        
