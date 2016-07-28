import json
import urllib2
import papersCache as pc
import re
import hashlib

def getDoi(doi):
  print "Getting DOI: "+doi
  check_doi = re.match(r'^http://dx\.doi\.org/', doi)
  if check_doi is not None:
    url = doi
    doi = re.sub(r'^http://dx\.doi\.org/', '', doi)
  else:
    url = 'http://dx.doi.org/'+doi
  request = urllib2.Request(url, headers={"accept": "application/json"})
  try:
    response = urllib2.urlopen(request)
    html_raw = response.read();
    json_data = json.loads(html_raw)
    #as doi's use '/' chars, we do an md5 of the doi as the filename
    filename = hashlib.md5(doi).hexdigest()
    pc.dumpJson(filename, json_data, filetype='raw/doi')
    return json_data
  except urllib2.HTTPError, e:
    print "DOI: "+doi+" error: "+str(e.code)
    return None

