import csv
import shutil

from config import config
from . import utils
from . import common_html as ch

###########################################################
# Build abstract word cloud
###########################################################
def build_abstract_word_cloud(data_from_count):

    print("\n###HTML - Abstract Word Cloud###")

    biggest_word_size = 100
    abstract_words = {}
    d3_word_list = "["
    n = 0

    with open(config.data_dir + "/abstract_raw.csv", 'rt') as f:

        reader = csv.reader(f)
        for row in reader:
            try:
                if row[0] != "":
                    abstract_words[row[0]] = int(row[1])
            except:
                pass

    # May not be any abstract text recieved (this can be domain specific)
    if len(abstract_words) == 0:
        return

    # Sort the words
    # abstract_words = sorted(abstract_words.items(), key=lambda x: x[1], reverse=True)
    # Sort by word first, then the count. This gives a consistent result as if the cut off is in a band
    # of words with the same counts it will randomly pick some. This pollutes commits to the web server
    abstract_words = sorted(dict(sorted(abstract_words.items())).items(), key=lambda x: x[1], reverse=True)

    # Only pick the top 200
    abstract_words = abstract_words[:200]

    # Grab the most frequent word
    max_count = abstract_words[0][1]

    # Scale the words so the biggest one is equal to biggest_word_size
    abstract_words = {k: round(v * biggest_word_size / max_count) for k, v in abstract_words}

    for this_word in sorted(abstract_words):
        if n > 0:
            d3_word_list += ","

        # word_list += '["' + row[0].replace("'","\'").replace('"','\"') + '",' + str(row[1]) +  ']'
        # word_list += '{"text":"' + row[0].replace("'", "\'").replace('"', '\"') + '","size":' + str(math.sqrt(int(row[1]))*1.5) + '}'
        # word_list += '{"text":"' + row[0].replace("'", "\'").replace('"', '\"') + '","size":' + str(row[1]) + '}'

        d3_word_list += '{"text":"' + this_word.replace("'", "\'").replace('"', '\"') + '","size":' + str(abstract_words[this_word]) + '}'
        n += 1

    d3_word_list += "];"
    list_file = open(config.html_dir + '/abstractwordcloud/list.js', 'w')
    list_file.write(" var word_list = " + d3_word_list)

    html_file = open(config.html_dir + '/abstractwordcloud/index.html', 'w')

    # Put html together for this page

    shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/abstractwordcloud/d3wordcloud.js')
    shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/abstractwordcloud/d3.layout.cloud.js')

    extra_head = '<script src="list.js"></script>'
    extra_head += '<script src="https://d3js.org/d3.v3.min.js"></script>'
    extra_head += '<script src="d3.layout.cloud.js"></script>'

    temp = ch.build_common_head("../", extra_head)
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Abstract Word Cloud</p>', "../")

    temp += '<h1 id="pagetitle">Abstract Word Cloud</h1>'

    temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
    temp += '</cloud>'

    temp += '<script src="d3wordcloud.js"></script>'
    temp += "<p>Data from " + utils.intWithCommas(data_from_count) + " publications. <span class='help_text'>(<a href='../help/index.html#missing_data'>What does this mean?</a>)</span></p>"

    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)


###########################################################
# Build keyword word cloud
###########################################################
def build_keyword_word_cloud(data_from_count):

    print("\n###HTML - Abstract Word Cloud###")

    biggest_word_size = 100
    words = {}
    d3_word_list = "["
    n = 0

    with open(config.data_dir + "/keywords_raw.csv", 'rt') as f:
        reader = csv.reader(f)
        for row in reader:
            try:
                if row[0] != "":
                    words[row[0]] = int(row[1])
            except:
                pass

    # Sort the words
    # words = sorted(words.items(), key=lambda x: x[1], reverse=True)
    # Sort by word first, then the count. This gives a consistent result as if the cut off is in a band
    # of words with the same counts it will randomly pick some. This pollutes commits to the web server
    words = sorted(dict(sorted(words.items())).items(), key=lambda x: x[1], reverse=True)

    # Only pick the top 200
    words = words[:200]

    # Grab the most frequent word
    max_count = words[0][1]

    # Scale the words so the biggest one is equal to biggest_word_size
    words = {k: round(v * biggest_word_size / max_count) for k, v in words}

    for this_word in words:
        if n > 0:
            d3_word_list += ","

        d3_word_list += '{"text":"' + this_word.replace("'", "\'").replace('"', '\"') + '","size":' + str(words[this_word]) + '}'
        n += 1

    d3_word_list += "];"
    list_file = open(config.html_dir + '/keyword_wordcloud/list.js', 'w')
    list_file.write(" var word_list = " + d3_word_list)

    html_file = open(config.html_dir + '/keyword_wordcloud/index.html', 'w')

    # # Put html together for this page

    shutil.copyfile(config.template_dir + '/d3wordcloud.js', config.html_dir + '/keyword_wordcloud/d3wordcloud.js')
    shutil.copyfile(config.template_dir + '/d3.layout.cloud.js', config.html_dir + '/keyword_wordcloud/d3.layout.cloud.js')

    extra_head = '<script src="list.js"></script>'
    extra_head += '<script src="https://d3js.org/d3.v3.min.js"></script>'
    extra_head += '<script src="d3.layout.cloud.js"></script>'

    temp = ch.build_common_head("../", extra_head)
    temp += ch.build_common_body('<p id="breadcrumbs"><a href="../index.html">Home</a> &gt; Keyword Word Cloud</p>', "../")

    temp += '<h1 id="pagetitle">Keyword Word Cloud</h1>'

    temp += '<cloud id="sourrounding_div" style="width:100%;height:500px">'
    temp += '</cloud>'

    temp += '<script src="d3wordcloud.js"></script>'
    temp += "<p>Data from " + utils.intWithCommas(data_from_count) + " publications. <span class='help_text'>(<a href='../help/index.html#missing_data'>What does this mean?</a>)</span></p>"

    html_file.write(temp)

    temp = ch.build_common_foot("../")
    html_file.write(temp)
