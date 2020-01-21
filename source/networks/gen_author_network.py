import json
# import unicodecsv
import os
import re
import sys


def loadCleaning():
    # Read in config file
    global pattern
    pattern = []
    global replacements
    replacements = []
    with open(os.path.join(os.getcwd(), '../config/institute_cleaning.csv'), 'rb') as csvfile:
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


if __name__ == '__main__':
    if __package__ is None:
        sys.path.append(os.getcwd())
        import config.config as config
        import analyse.genLinks as gl
    else:
        from ..analyse import genLinks as gl
        from ..config import config as config
else:
    import analyse.genLinks as gl
    import config.config as config


def output_network():
    # setup a new instance of dataNetwork, specifying source of data
    dn = gl.dataNetwork(os.path.join(config.cache_dir, 'processed/merged'))
    loadCleaning()
    dn.target_path = '$.merged.author.[*]'
    dn.node_name = ['$.family', '$.given']
    dn.additional_data_node.append({'path': '$.affiliation.[0].name', 'context': 'node', 'name': 'affiliation', 'clean': cleanInstitution})
    dn.genNetwork()
    output = dn.processOutput()
    f = open(os.path.join(config.cache_dir, 'processed/authorlinks.cleaned.json'), 'w')
    json.dump(output, f, indent=2)
    f.close()


if __name__ == '__main__':
    # Lets figure out some paths that everything is relative to
    # global root_dir
    path_to_papers_py = os.path.abspath(sys.argv[0])
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(path_to_papers_py)))
    print('Root directory = ' + root_dir)

    # Get all the config - these will be a global vars available like config.varname
    config.build_config_variables(root_dir)

    output_network()
