import re
import json
import sys
import os
import config.config as config


# dump dict data to filename (cache[/filetype]/filename
# if process is true adds '[' & ']' to start and end otherwise json will not parse correctly
def dumpJson(filename, data, filetype='', process=False):
    try:
        location = '/'.join(filter(None, [config.cache_dir, filetype, filename]))
        f = open(location, 'w')
        if process:
            f.write('[')
            for datum in data:
                json.dump(datum, f, indent=4)
                f.write(',')

            f.seek(-1, os.SEEK_END)
            f.truncate()
            f.write(']')
        else:
            json.dump(data, f, indent=4)
        f.close()
        return location
    except:
        print("(papersCache.dumpJson) Unexpected error:", sys.exc_info()[1])
        pass


# dump dict data to filename (cache[/filetype]/filename
def dumpFile(filename, data, filetype=''):
    try:
        location = '/'.join(filter(None, [config.cache_dir, filetype, filename]))
        print("location:" + str(type(location)))
        print(location)
        print("data:" + str(type(data)))
        f = open(location, 'w')
        f.write(data)
        f.close()
        return location
    except:
        print("(papersCache.dumpFile) Unexpected error:", sys.exc_info()[1])
        pass


# get a list of all filenames in directory cache[/filetype]
def getCacheList(filetype=''):
    cachefiles = []
    filetype = re.sub('[.]{2,}', '', filetype)
    for root, dirs, files in os.walk(config.cache_dir+filetype):
        for name in files:
            cachefiles.append(name)
    return cachefiles


# get all [filenames] from cache[/filetype]
# processes to output_type (def. json)
# returns a dict of [filename]->[filedata]
def getCacheData(output_type='json', filetype='', filenames=[]):
    cache_files = []
    location = '/'.join(filter(None, [config.cache_dir, filetype]))
    for root, dirs, files in os.walk(config.cache_dir+'/'+filetype):
        for name in files:
            if (len(filenames) == 0):
                cache_files.append(name)
            else:
                if name in filenames:
                    cache_files.append(name)

    papers = {}

    for cache_filename in cache_files:
        location = '/'.join(filter(None, [config.cache_dir, filetype, cache_filename]))
        cache_file = open(location, 'r')
        # cache_file_str = ''.join(cache_file.read().split())
        # print cache_file_str
        # paper = json.loads(cache_file_str)
        if (output_type == 'json'):
            papers[cache_filename] = json.load(cache_file)
            cache_file.close()
        else:
            raise ValueError('(papersCache.getCacheData) Error, unrecognised output_type')

    return papers
