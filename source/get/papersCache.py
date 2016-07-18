import re
import hashlib
import json
import sys
import os
import time

#dump dict data to filename (cache[/filetype]/filename
#if process is true adds '[' & ']' to start and end otherwise json will not parse correctly
def dumpJson(filename, data, filetype='', process=False):
  try:
    location = '/'.join(filter(None, ['../../cache', filetype, filename]))
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

#dump dict data to filename (cache[/filetype]/filename
def dumpFile(filename, data, filetype=''):
  try:
    location = '/'.join(filter(None, ['../../cache', filetype, filename]))
    f = open(location, 'wa')
    f.write(data)
    f.close()
    return f
  except:
    print "(papersCache.dumpJson) Unexpected error:", sys.exc_info()[1]
    pass

#get a list of all filenames in directory cache[/filetype]
def getCacheList(filetype = ''):
  cachefiles = []
  filetype = re.sub('[\.]{2,}', '', filetype)
  for root, dirs, files in os.walk('../../cache/'+filetype):
    for name in files:
      cachefiles.append(name)
  return cachefiles

#get all [filenames] from cache[/filetype]
#processes to output_type (def. json)
#returns a dict of [filename]->[filedata]
def getCacheData(output_type='json', filetype='', filenames=[]):
  cache_files = []
  location = '/'.join(filter(None, ['../../cache', filetype]))
  for root, dirs, files in os.walk('../../cache/'+filetype):
    for name in files:
      if (len(filenames) == 0):
        if re.search('^WHAT GOES HERE?$', name) != None:
          cache_files.append(name)
      else:
        if name in filenames:
          cache_files.append(name)
   
  papers = {}

  for cache_filename in cache_files:
    location = '/'.join(filter(None, ['../../cache', filetype, cache_filename]))
    cache_file = open(location, 'r')
    #cache_file_str = ''.join(cache_file.read().split())
    #print cache_file_str
    #paper = json.loads(cache_file_str)
    if (output_type == 'json'):
      papers[cache_filename] = json.load(cache_file)
    else:
      raise ValueError('(papersCache.getCacheData) Error, unrecognised output_type')

  return papers    
