import re
import hashlib
import json
import sys
import os
import time

def dumpJson(filename, data, filetype='', process=True):
  try:
    location = '/'.join(filter(None, ['./cache', filetype, filename]))
    f = open(location, 'wa')
    if process:
      f.write('[')
      for datum in data:
        json.dump(datum, f, indent=2)
        f.write(',')

      f.seek(-1, os.SEEK_END)
      f.truncate()
      f.write(']')
    else:
      json.dump(data, f, indent=2)
    f.close()
    return f
  except:
    print "(papersCache.dumpJson) Unexpected error:", sys.exc_info()[1]
    pass

def getCacheList(filetype = ''):
  cache_files = []
  filetype = re.sub('[\.]{2,}', '', filetype)
  for root, dirs, files in os.walk('./cache/'+filetype):
    for name in files:
      cachefiles.append(name)
  return cachefiles

def getCacheData(output_type='json', filetype='', filenames=[]):
  cache_files = []
  for root, dirs, files in os.walk('./cache/'+filetype):
    for name in files:
      if (len(filenames) == 0):
        if re.search('^WHAT GOES HERE?$', name) != None:
          cache_files.append(name)
      else:
        if name in filenames:
          cache_files.append(name)
   
  papers = {}

  for cache_filename in cache_files:
    cache_file = open('./cache/'+filetype+'/'+cache_filename, 'r')
    #cache_file_str = ''.join(cache_file.read().split())
    #print cache_file_str
    #paper = json.loads(cache_file_str)
    if (output_type == 'json'):
      papers[cache_filename] = json.load(cache_file)
    else:
      raise ValueError('(papersCache.getCacheData) Error, unrecognised output_type')

  return papers    
