#authors (first, simple, with affiliation)

import analyse.genLinks as gl
import json
import unicodecsv
import re

def loadCleaning():
  # Read in config file
  global pattern
  pattern = []
  global replacements
  replacements = []
  with open('/home/lishy/repos/papers.origin/config/institute_cleaning.csv', 'rb') as csvfile:
    f = unicodecsv.reader(csvfile,  encoding='utf-8')  # Handle extra unicode characters
    for row in f:
      try:
        # Check it is not a comment string first.
        if re.match('#', row[0]):
          continue

        # Check for blank lines
        if row[0] == '':
          continue

        # If there is a second element in this row then carry on
        pattern.append(row[0])
        replacements.append(row[1])
      except:
        pass

def cleanInstitution(value):
  # Cycle through institute checking the whole substitution list.
  # Stop when the first one matches.
  for y in range(0, len(pattern)):
    temp = re.search(pattern[y], value, re.IGNORECASE)
    if temp > 0:
      value = replacements[y]
      break

  return value

dn = gl.dataNetwork('/home/lishy/repos/papers.origin/cache/alspac/processed/merged')
loadCleaning()
dn.target_path = '$.author.[0]'
dn.node_name = ['$.family', '$.given']
dn.additional_data_node.append({'path': '$.affiliation.[0].name', 'context': 'node', 'name':'affiliation', 'clean': cleanInstitution})
dn.genNetwork()
output = dn.processOutput()
f = open('/home/lishy/repos/papers.origin/cache/alspac/processed/firstauthorlinks.cleaned.json', 'w')
json.dump(output, f, indent=2)
f.close()

