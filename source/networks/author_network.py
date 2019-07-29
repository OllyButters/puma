import config.config as config
import gen_author_network
import gen_html
import os


def build_network():

    # first generate the network datafile
    gen_author_network.output_network()

    # this dumps the output to processed/authorlinks.cleaned.json

    datafile_name = 'processed/authorlinks.cleaned.json'
    title = 'Author network'
    # script name refers to script to use from html/templates
    script_name = 'network.authors.js'

    # outputdir is where to dump generated html and javascript files
    # root path is config.html_dir
    outputdir = 'authornetwork'

    # generate the html (uses a template in html/templates
    network_html = gen_html.build_network_page(network_title=title, network_datafile=os.path.split(datafile_name)[1], script_name=script_name)

    # copy across the css and js files required
    gen_html.copy_network_files(outputdir=outputdir, network_datafile=datafile_name, script_name=script_name)

    with open(os.path.join(config.html_dir, outputdir, os.path.split(datafile_name)[1] + '.html'), 'w') as outputfile:
        outputfile.write(network_html)
