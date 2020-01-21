import json
#import unicodecsv
import re
import os
import getopt
import sys
from pprint import pprint
from string import Template
import shutil


def copy_network_files(network_datafile, script_name, outputdir):
  # Copy CSS files
  shutil.copyfile(config.template_dir + '/style_main.css', config.html_dir + '/css/style_main.css')

  #copy js file
  shutil.copyfile(os.path.join(config.template_dir, script_name), os.path.join(config.html_dir, outputdir, script_name))

    #copy images
  imagefiles = (
    'dot-highlight.png',
    'vline3.png',
    'point-40-1px.png',
    'dot4.png',
    'dot.png',
  )
  for imagefile in imagefiles:
    shutil.copyfile(os.path.join(config.template_dir, imagefile), os.path.join(config.html_dir, outputdir, imagefile))



  #copy datafile
  shutil.copyfile(os.path.join(config.cache_dir, network_datafile), os.path.join(config.html_dir, outputdir, os.path.split(network_datafile)[1]))

def build_network_page(network_title, network_datafile, script_name):
  #load the html template file
  with open(os.path.join(config.template_dir, 'network.template.pixi.html'), 'r') as template_file:
    template_text = template_file.read()

  template = Template(template_text)

  #add in html common body
  html_text = template.safe_substitute(common_body = build_html.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; ' + network_title + '</p>', "../", ""), common_foot = build_html.build_common_foot())

  template = Template(html_text)
  html_text = template.safe_substitute(title = network_title, datafile = network_datafile, script_name = script_name)

  return html_text


if __name__ == '__main__':
  if __package__ is None:
    sys.path.append(os.getcwd())
    import config.config as config
    import html.build_htmlv2 as build_html
  else:
    from ..config import config as config
    from ..html import build_htmlv2 as build_html
else:
  import web_pages.build_htmlv2 as build_html
  import config.config as config


if __name__ == '__main__':
  # Lets figure out some paths that everything is relative to
  # global root_dir
  path_to_papers_py = os.path.abspath(sys.argv[0])
  root_dir = os.path.dirname(os.path.dirname(os.path.dirname(path_to_papers_py)))
  print('Root directory = ' + root_dir)

  # get the args relevant for script inputs
  try:
    opts, args = getopt.getopt(sys.argv[1:], "o:d:t:h:", ["outputdir=", "datafile=", "net_type=", "title="])
  except Exception as e:
    pprint(str(e))
    sys.exit(2)
  # only pass config arg to config.ini
  sys.argv = [sys.argv[0]] + [[v, sys.argv[i+1]] for i,v in enumerate(sys.argv) if v == 'config' or v == 'c']


  # Get all the config - these will be a global vars available like config.varname
  config.build_config_variables(root_dir)

  # sort out other options
  for opt, arg in opts:
    if opt in ('o', '--outputdir'):
      outputdir = arg
    elif opt in ('d', '--datafile'):
      datafile = arg
    elif opt in ('t', '--net_type'):
      net_type = arg
    elif opt in ('h', '--title'):
      title = arg

  #add in script name
  if net_type == 'simple_nodes':
    script_name = 'network.simplenodes.js'
  elif net_type == 'authors':
    script_name = 'network.authors.js'
  else:
    script_name = 'network.simplenodes.js'

  datafile_name = os.path.split(datafile)[1]

  network_html = build_network_page(network_title = title, network_datafile = datafile_name, script_name = script_name)

  copy_files(outputdir = outputdir, network_datafile = datafile, script_name = script_name)

  with open(os.path.join(config.html_dir, outputdir, datafile_name + '.html'), 'w') as outputfile:
    outputfile.write(network_html)
