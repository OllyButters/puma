import re
import json
import papersCache as pc
import jsonpath_rw as jsonp
import pprint
import copy
import logging

class Merge():
  def __init__(self):
    self.src = None
    self.dest = None
    self.mapping = {}
    self.output = []
    self.log = logging.getLogger('mergeLog')
    self.log.setLevel(logging.INFO)
    fh = logging.FileHandler(filename='../../logs/mergeLog', mode="w")
    self.log_fh = fh
    fh.setLevel(logging.INFO)
    fh_format = logging.Formatter('%(message)s')
    fh.setFormatter(fh_format)
    self.log.addHandler(fh)
    self.log.info("Merge INFO - class instantiated")

  def __del__(self):
    self.log.removeHandler(self.log_fh)
    self.log_fh.close()

  #iterate over all fields in src_data (nested dict/list) and append to dest if it doesn't currently exist
  #recursively called for sub- lists/dicts
  #todo needs rationalising to remove searches for matches. data can just be added by matching dest path or src path if not exist
  def iterFields(self, src_data, path_mapping, full_src_path, dest, src_path = None, dest_parent_path = None):

    #if src_path is None, set it to full_src_path
    if src_path is None:
      src_path = full_src_path

    path_mapping_data = self.findDestPath(path_mapping, full_src_path)
    dest_path = path_mapping_data['path']
    path_subfields = path_mapping_data['sub_paths']
    dest_path_last_element = dest_path.split('.')[-1]

    #assign dest_parent_path to new_dest_parent_path
    new_dest_parent_path = dest_path

    #so now we know what (if any) mappings occur for this path, we either need to point to the right dict/list in
    #dest or add as new and point to the location in dest
    #first we try to find it
    if dest_path is not None and dest_path != '':
      #trim off the destination_parent_path
      if dest_parent_path is not None and dest_parent_path != '' and dest_path.startswith(dest_parent_path):
        dest_path = '$'+dest_path[len(dest_parent_path):]

      matches = self.getPath(dest, dest_path)

      src_loc_type = re.match(r'\[([0-9]+)\]', src_path[-1])

      dest_loc = -1
      if src_loc_type is not None:
        dest_loc = int(src_loc_type.group(1))
      
      if len(matches) > 0:
        #we've found some matches
        #self.log.info('Some matches found in dest: '+str(matches))
        #what type is the match?
        if isinstance(matches[dest_loc].value, dict):
          #if it's a dict, we can just set it to value
          dest = matches[dest_loc].value
          if isinstance(matches[dest_loc].path, jsonp.Root):
            dest_field = None
          else:
            dest_field = str(matches[dest_loc].path)
        elif isinstance(matches[dest_loc].value, list):
          dest = matches[dest_loc].value
          if len(dest) == 0:
            dest.append('')
          dest_field = len(matches[dest_loc].value) - 1
        else:
          #it's a str or int (ignore tuples/sets/etc)
          dest = matches[dest_loc].context.value
          dest_field = matches[dest_loc].path
      else:
        #if no matches found, test to see if in self.dest as opposed to dest. if not there either, insert new location and set as dest
        #note that new_dest_parent_path is dest_path before trimming of dest_parent_path
        #todo rationalise this into single function
        matches = self.getPath(self.dest, new_dest_parent_path)

        if len(matches) > 0:
          #we've found some matches
          #self.log.info('Some matches found in self.dest: '+str(matches))
          #what type is the match?
          if isinstance(matches[-1].value, dict):
            #if it's a dict, we can just set it to value
            dest = matches[-1].value
            if isinstance(matches[-1].path, jsonp.Root):
              dest_field = None
            else:
              dest_field = str(matches[-1].path)
          elif isinstance(matches[-1].value, list):
            dest = matches[-1].value
            if len(dest) == 0:
              dest.append('')
            dest_field = len(matches[-1].value) - 1
          else:
            #it's a str or int (ignore tuples/sets/etc)
            dest = matches[-1].context.value
            dest_field = matches[-1].path
        else:
          #no matches found, so go on to add full_src_path to self.dest
          dest_data = self.insertNewDestLocation(full_src_path)
          dest = dest_data['dest']
          dest_field = dest_data['dest_field']
          dest_path = '.'.join(src_path)
    elif dest_path == '':
      #if the destination path doesn't exist, insert new location as above and set as dest
      dest_data = self.insertNewDestLocation(full_src_path)
      dest = dest_data['dest']
      dest_field = dest_data['dest_field']
      dest_path = '.'.join(src_path)

    #what type is data?
    #if it's a dict/list, we want to go through the items by key/index and check type
    #if it's not the above (i.e. str or int (ignore sets/tuples)) then we need to find out where the PATH maps to and then set the DEST value to data
    if isinstance(src_data, list) or isinstance(src_data, dict):
      #ok, so it's a list or dict

      if isinstance(dest, list):
        if isinstance(src_data, list):
          while len(dest) < len(src_data):
            dest.append(copy.deepcopy(dest[0]))
        #elif isinstance(src_data, dict):
        #  while len(dest) < len(src_data.keys()):
        #    dest.append(copy.deepcopy(dest[0]))
 
      if isinstance(src_data, list):
        for key, item in enumerate(src_data):
          #the sub_src_path is effectively an id for this particular entry in data (source data)
          sub_src_path = full_src_path + ('['+str(key)+']',)
          src_path = ['$', '['+str(key)+']']
          #now we call iterfields again, setting dest as the found match, or created entry and the src_data, src_path and dest_path as required
          #if dest is a list, we need to copy the last sub-element if present and we're not on the last element of src_data
          #if isinstance(dest, list):
          #  if len(dest) > 0:
          #    #if key < len(src_data) - 1:
          #    dest_field = len(dest) - 1
          #    dest.append(copy.deepcopy(dest[0]))
          #      #print "dest is list"
          #    #else:
          #    if key == len(src_data) - 1:
          #      dest.pop(0)
          #    #  if len(dest) > 1:
          #    #    dest.pop(0)
          #    #  dest_field = len(dest) - 1
          self.iterFields(item, path_subfields, sub_src_path, dest, src_path, new_dest_parent_path)
          #if isinstance(dest, list):
          #  dest.pop(0)
      elif isinstance(src_data, dict):
        for fieldname, item in src_data.items():
          #again, the sub_src_path is effectively an id for this particular entry in data (source data)
          sub_src_path = full_src_path + (fieldname,)
          src_path = ['$', fieldname]
          #now we call iterfields again, setting dest as the found match, or created entry and the src_data, src_path and dest_path as required
          #self.log.info('--end iterFields-----------------------------------------')
          #self.iterFields(item, path_subfields, sub_src_path, dest, src_path, dest_path)
          self.iterFields(item, path_subfields, sub_src_path, dest, src_path, new_dest_parent_path)
    else:
      #data is str, int, etc (n.b. this ignores sets, tuples as they don't currently exist in data, but may be a future issue)
      #we therefore assign value to dest[key]
      if not(isinstance(src_data, str) and src_data == '') and src_data != None:
        if isinstance(dest, list):
          if isinstance(dest_field, jsonp.Index):
            dest_field = re.sub('[\[\]]', '', str(dest_field))
          dest[int(dest_field)] = src_data
          #try:
          #  dest[int(dest_field)] = src_data
          #except IndexError:
          #  for i in range(len(dest), int(dest_field) + 1):
          #    dest.append({})
          #  dest[int(dest_field)] = src_data
        elif isinstance(dest, dict):
          dest[str(dest_field)] = src_data
        else:
          dest = src_data

  #mapSrc calls iterfields intially, setting the correct values to begin the mapping
  def mapSrc(self):
    full_path = ('$',)
    self.iterFields(self.src, self.mapping, full_path, self.dest)

  #try to find path in path_mapping, return a dict of the path (if found) and any sub-paths
  def findDestPath(self, path_mapping, path):
    json_path = '.'.join(path)
    json_path = re.sub('\[[0-9]+\]', '[*]', json_path)
    
    lookup_path = ""
    sub_paths = {}
    
    if json_path in path_mapping.keys():
      lookup_path = path_mapping[json_path]
      for key, path in path_mapping.items():
        if key.startswith(json_path) and key != json_path:
          sub_paths[key] = path
    else:
      sub_paths = path_mapping

    mapping_output = {
      'path': lookup_path,
      'sub_paths': sub_paths
    }

    return mapping_output

  def getPath(self, target_data, path):
    #parse the path as json
    json_path = jsonp.parse(path)
    #create a list of matches
    matches = [match for match in json_path.find(target_data)]
    #return the list of matches (DatumInContext objects, key properties: value, context (another DatumInContext object)) 
    #when we assign the value it must be in the jsonpath_rw match object.context.value.[final part of path] (DatumInContext.context.value.[pathfield])
    return matches

  #insert path into destination
  def insertNewDestLocation(self, path, dest = None):
    if dest is None:
      dest = self.dest
    current_path = dest
    #path includes the '$' character as the first element, so remove this
    path = path[1:]
    if len(path) == 0:
      return {'dest': self.dest, 'dest_field': None}
    new_path = []

    #we need to go through path and construct the whole new path
    for loc in path:
      #if loc matches the list format pattern, we need to split it into two elements: the list ([]) and the index (n)
      loc_type = re.match(r'\[([0-9]+)\]', loc)

      if loc_type is not None:
        new_path.append([])
        new_path.append(loc_type.group(1))
      else:
        new_path.append(loc)

    #we now have a new path constructed of strings and list elements (e.g. ['foo', 'bar', [], '0', 'car'])
    dest_field = new_path[-1]
    #if the last elements are a list and then an index, we need to drop the index from the end otherwise we end up pointing at the wrong place
    if len(new_path) > 1 and isinstance(new_path[-2], list):
      new_path = new_path[:-1]

    i = 0

    while i < len(new_path) - 1:
      loc = new_path[i]
      #print "loc: "+str(loc)
      if i + 1 > len(new_path) - 1:
        next_loc = None
      else:
        next_loc = new_path[i+1]

      if isinstance(next_loc, list):
        next_element = []
      else:
        next_element = {}

      if loc is not None:
        if isinstance(current_path, list):
          #we don't need to test for loc being a list here as it must always be an int (type str) (see splitting path above)
          if len(current_path) <= int(next_loc):
            current_path.append(next_element)
            current_path = current_path[-1]
            #we need to skip adding the next element (always an index following a list element) otherwise it will get added as a dict key
            i += 1
          else:
            current_path = current_path[int(next_loc)]
            #we need to skip adding the next element (always an index following a list element) otherwise it will get added as a dict key
            i += 1
        else:
          if loc not in current_path:
            current_path[loc] = next_element
            current_path = current_path[loc]
          else:
            current_path = current_path[loc]
      i += 1

    if isinstance(current_path, list):
      if len(current_path) <= int(dest_field):
        current_path.append({})

    return {'dest': current_path, 'dest_field': dest_field}

  #nicked from https://www.xormedia.com/recursively-merge-dictionaries-in-python/
  def dict_merge(self, a, b):
    '''recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and bhave a key who's value is a dict then dict_merge is called
    on both values and the result stored in the returned dictionary.'''
    #note that this is adjusted to use lists as well as dicts. the type of object in a will override that of b.
    if not (isinstance(b, dict) or isinstance(b, list)):
      return b
    result = copy.deepcopy(a)
    if isinstance(result, dict):
      res_type = 'dict'
    elif isinstance(result, list):
      res_type = 'list'
    else:
      res_type = 'other'

    if isinstance(b, dict):
      for k, v in b.iteritems():
        check = False
        if res_type == 'dict':
          check = k in result
        else:
          check = True if len(result) > k else False
          if not check:
            result.append('')

        if check and (isinstance(result[k], dict) or isinstance(result[k], list)):
          result[k] = self.dict_merge(result[k], v)
        else:
          if v != "":
            result[k] = copy.deepcopy(v)
    else:
      #b is a list
      for k, v in enumerate(b):
        check = False
        if res_type == 'dict':
          check = k in result
        else:
          check = True if len(result) > k else False
          if not check:
            result.append('')

        if check and (isinstance(result[k], dict) or isinstance(result[k], list)):
          result[k] = self.dict_merge(result[k], v)
        else:
          if v != "":
            result[k] = copy.deepcopy(v)

    return result

