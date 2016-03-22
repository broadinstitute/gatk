import requests
import sys
import time

ip="146.148.89.97"
url = 'http://'+ip+':8000/api/workflows/v1'
wdl = 'hello.wdl'
files = {'wdlSource': open(wdl, 'rb')}

def launchHelloWdl():
    r = requests.post(url, files = files)
    r.raise_for_status()
    return r.json()['id']

def checkCompletion(id):
    r = requests.get("http://"+ip+":8000/api/workflows/v1/"+id+"/status")
    r.raise_for_status()
    return r.json()["status"]

def main(args):
    id = launchHelloWdl()
    done = False
    while not done:
        result = checkCompletion(id)
        print result
        time.sleep(1)
        if result == "Succeeded":
            done = True


    print "done"

if __name__ == "__main__":
    main(sys.argv[1:])

