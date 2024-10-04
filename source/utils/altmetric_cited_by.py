"""Use DOIs to get counts of citations from the Altmetric API. These are citations like tweets, blog posts, etc."""

import argparse
import csv
import json
import os
import requests
import sys

CACHE_ROOT = 'cache/'
BASE_URL = 'https://api.altmetric.com/v1/'

def get_altmetric_data(doi_list):
    """Get the Altmetric data for a list of DOIs from the Altmetric API and save it as JSON files."""

    for doi in doi_list:
        print('Getting data for doi:', doi)

        url = BASE_URL + 'doi/' + doi

        try:
            response = requests.get(url, timeout=10)
            print("Response code: " + str(response.status_code))

            if response.status_code == 200:

                try:
                    data = response.json()
                    output_filename = doi.replace('/', '_') + '.json'
                    output_filename = CACHE_DIR + output_filename
                    with open(output_filename, 'w', encoding="utf-8") as f:
                        f.write(json.dumps(data, indent=4))

                except Exception as e:
                    print('Error with doi:', doi)
                    print(e)
                    continue

        except response.exceptions.RequestException as e:
            print('Error with doi:', doi)
            print(e)


def parse_altmetric_data():
    """Parse the Altmetric data from the JSON files and write it to a CSV file."""
    input_files = os.listdir(CACHE_DIR)
    output_file = 'output_' + project + '.csv'

    with open(output_file, 'w', newline='', encoding="utf-8") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Title',
                            'Cited by Facebook Walls', 
                            'Cited by Blogs', 
                            'Cited by Mainstream Media', 
                            'Cited by Redditors', 
                            'Cited by Forums and Stack Exchange',
                            'Cited by Twitter accounts',  
                            'Cited by Wikipedia', 
                            'Cited by Policies', 
                            'Cited by Patents', 
                            'Cited by YouTube channels',
                            'Cited by Total number of accounts',
                            'Cited by Total number of posts'])

        for input_file in input_files:
            with open(CACHE_DIR + input_file, 'r', encoding="utf-8") as f:
                data = json.load(f)
                print(data)
                try:
                    title = data['title']
                except KeyError:
                    title = 'No title'

                try:
                    cited_by_fbwalls_count = data['cited_by_fbwalls_count']
                except KeyError:
                    cited_by_fbwalls_count = 0

                try:
                    cited_by_feeds_count = data['cited_by_feeds_count']
                except KeyError:
                    cited_by_feeds_count = 0

                try:
                    cited_by_msm_count = data['cited_by_msm_count']
                except KeyError:
                    cited_by_msm_count = 0

                try:
                    cited_by_rdts_count = data['cited_by_rdts_count']
                except KeyError:
                    cited_by_rdts_count = 0

                try:
                    cited_by_qna_count = data['cited_by_qna_count']
                except KeyError:
                    cited_by_qna_count = 0

                try:
                    cited_by_tweeters_count = data['cited_by_tweeters_count']
                except KeyError:
                    cited_by_tweeters_count = 0

                try:
                    cited_by_wikipedia_count = data['cited_by_wikipedia_count']
                except KeyError:
                    cited_by_wikipedia_count = 0

                try:
                    cited_by_policies_count = data['cited_by_policies_count']
                except KeyError:
                    cited_by_policies_count = 0

                try:
                    cited_by_patents_count = data['cited_by_patents_count']
                except KeyError:
                    cited_by_patents_count = 0

                try:
                    cited_by_videos_count = data['cited_by_videos_count']
                except KeyError:
                    cited_by_videos_count = 0

                try:
                    cited_by_accounts_count = data['cited_by_accounts_count']
                except KeyError:
                    cited_by_accounts_count = 0

                try:
                    cited_by_posts_count = data['cited_by_posts_count']
                except KeyError:
                    cited_by_posts_count = 0

                csvwriter.writerow([title,
                                    cited_by_fbwalls_count,
                                    cited_by_feeds_count,
                                    cited_by_msm_count,
                                    cited_by_rdts_count,
                                    cited_by_qna_count,
                                    cited_by_tweeters_count,
                                    cited_by_wikipedia_count,
                                    cited_by_policies_count,
                                    cited_by_patents_count,
                                    cited_by_videos_count,
                                    cited_by_accounts_count,
                                    cited_by_posts_count])



# Get the input file and project short name from the command line
parser = argparse.ArgumentParser()
parser.add_argument("--input", help="CSV file with DOIs")
parser.add_argument("--project", help="Project short name")
args = parser.parse_args()

if args.input:
    print('CSV DOI file set to: '+args.input)
    input_file_name = args.input
else:
    print('No input file set! Use --input to set the input file.')
    sys.exit()

if not os.path.isfile(input_file_name):
    print('Input file does not exist!')
    sys.exit()

if args.project:
    print('Project short name set to: '+args.project)
    project = args.project
else:
    print('No project short name set! Use --project to set the project short name.')
    sys.exit()

# Make sure the cache directory exists
CACHE_DIR = CACHE_ROOT + "/" + project + '/'
os.makedirs(CACHE_DIR, exist_ok=True)


doi_list = []

with open(input_file_name, 'r', encoding="utf-8") as f:
    reader = csv.reader(f)
    header = next(reader)
    for row in reader:
        this_doi = row[8]
        doi_list.append(this_doi)

print(doi_list)

get_altmetric_data(doi_list)
parse_altmetric_data()
