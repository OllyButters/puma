# Publications Metadata Augmentation (PUMA) pipeline

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3971102.svg)](https://doi.org/10.5281/zenodo.3971102)

PUMA is a python pipeline which aggregates third party metadata about journal articles to enable bibliometric research and easy exploration and searching of a collection of articles.

The technical documentation with screenshots is at: [https://github.com/OllyButters/puma/wiki](https://github.com/OllyButters/puma/wiki)

A demo with some of the functionality is at: [https://ollybutters.github.io/puma_web](https://ollybutters.github.io/puma_web)

## Running the code

The wiki (<https://github.com/OllyButters/puma/wiki>) has the best instructions, but in overview:

### Linux

- Assuming a linux environment with python 3 installed.
- Add all the required python libraries with the requirements.txt file in the source folder:
  
    `pip3 install -r requirements.txt`

- Edit the config.ini file in the config folder so it has your API keys etc.
- Run the pipeline from within the source folder with a command like:
  
    `./papers.py`

- Or you can run it from within an IDE.
- Once it has finished, have a look at the web pages in the html folder. Pay particular attention to the coverage_report.html file.

### Windows

- Assuming a windows environment (I've only tried Windows 10) with python 3 installed (at least python 3.6 as utf-8 support is a little flakey before that).
- Add all the required python libraries with the requirements.txt file in the source folder by running a command like the following with the command prompt:
  
    `py -3 -m pip install --user -r requirements.txt`
- Edit the config.ini file in the config folder so it has your API keys etc.
- Run the pipeline with a command like:

    `python papers.py`

- Or you can run it from within an IDE.
- Once it has finished, have a look at the web pages in the html folder. Pay particular attention to the coverage_report.html file.

### Azure pipelines

The example instances of the output at [https://ollybutters.github.io/puma_web](https://ollybutters.github.io/puma_web) are automatically built and deployed with the `azure-pipelines.yml` file in the root of the web repo at <https://github.com/OllyButters/puma_web/>.

## PUMA project publications

- Butters et al., 2019. DOI: [10.12688/wellcomeopenres.14986.1](https://dx.doi.org/10.12688%2Fwellcomeopenres.14986.1)
- Butters et al., 2021. DOI: [10.12688/f1000research.25484.2](https://doi.org/10.12688/f1000research.25484.2)

## Project Leads

- Dr Becca Wilson, University of Liverpool Co-I & project founder
- Dr Olly Butters, University of Liverpool Co-I

## Contributors

- Hugh Garner, University of Newcastle
- Tom Burton, University of Bristol
- Amran Ismail, University of Bristol
- Sue Thompson, Newcastle University

## Funding

Puma has been running for a couple of years and has been funded by a few different sources, of particular note are:

- Nuffield Foundation Summer Research placement
- D2K, University of Bristol
- [CLOSER](https://closer.ac.uk)
- [D2K, Newcastle University](https://research.ncl.ac.uk/d2k/)
- UKRI Innovation Fellowship with HDR UK
- [NIHR ARC NWC](https://arc-nwc.nihr.ac.uk/)
