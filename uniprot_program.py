#!/usr/bin/env python3

# standard library modules
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep

BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/entry/InterPro/protein/UniProt/?page_size=100"

def output_list():
  #disable SSL verification to avoid config issues
  context = ssl._create_unverified_context()

  next_ = BASE_URL
  last_page = False

  
  attempts = 0
  while next_:
    try:
      req = request.Request(next_, headers={"Accept": "application/json"})
      res = request.urlopen(req, context=context)
      # If the API times out due a long running query
      if res.status == 408:
        # wait just over a minute
        sleep(61)
        # then continue this loop with the same URL
        continue
      elif res.status == 204:
        #no data so leave loop
        break
      payload = json.loads(res.read().decode())
      next_ = payload["next_"]
      attempts = 0
      if not next_:
        last_page = True
    except HTTPError as e:
      if e.code == 408:
        sleep(61)
        continue
      else:
        # If there is a different HTTP error, it wil re-try 3 times before failing
        if attempts < 3:
          attempts += 1
          sleep(61)
          continue
        else:
          sys.stderr.write("LAST URL: " + next_)
          raise e

    for i, item in enumerate(payload["results"]):
      
        sys.stdout.write(item["metadata"]["accession"] + "\n")
      
      # Don't overload the server, give it time before asking for more
    if next_:
      sleep(1)

if __name__ == "__main__":
  output_list()

