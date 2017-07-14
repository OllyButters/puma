#example of associations across papers between mesh
#headings and subjects 

import json
import unicodecsv
import re
import os
import sys
from pprint import pprint

if __name__ == '__main__':
  if __package__ is None:
    print 'getcwd: '+os.getcwd()
    sys.path.append(os.getcwd())
    import config.config as config
    import analyse.genLinks as gl
  else:
    from ..analyse import genLinks as gl
    from ..config import config as config

# Lets figure out some paths that everything is relative to
# global root_dir
path_to_papers_py = os.path.abspath(sys.argv[0])
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(path_to_papers_py)))
print 'Root directory = ' + root_dir

# Get all the config - these will be a global vars available like config.varname
config.build_config_variables(root_dir)

#setup a new instance of dataNetwork, specifying source of data
dn = gl.dataNetwork(os.path.join(config.cache_dir, 'processed/merged'))
#get the links between the two sets of nodes
dn.target_path = '$.merged.subject.[*]'
dn.assoc_node_target_path = '$.merged.MedlineCitation.MeshHeadingList.[*].DescriptorName'
dn.search_across = 'all_data'
dn.genNetwork()
output = dn.processOutput()
f = open(os.path.join(config.cache_dir, 'processed/subject.mesh.links.json'), 'w')
json.dump(output, f, indent=2)
f.close()

#now run standard links across papers for nodes and assoc_nodes
dn = gl.dataNetwork(os.path.join(config.cache_dir, 'processed/merged'))
dn.target_path = '$.merged.subject.[*]'
dn.search_across = 'paper'
dn.genNetwork()
output = dn.processOutput()
f = open(os.path.join(config.cache_dir, 'processed/subject.links.json'), 'w')
json.dump(output, f, indent=2)
f.close()

dn = gl.dataNetwork(os.path.join(config.cache_dir, 'processed/merged'))
dn.target_path = '$.merged.MedlineCitation.MeshHeadingList.[*].DescriptorName'
dn.search_across = 'paper'
dn.genNetwork()
output = dn.processOutput()
f = open(os.path.join(config.cache_dir, 'processed/mesh.links.json'), 'w')
json.dump(output, f, indent=2)
f.close()
