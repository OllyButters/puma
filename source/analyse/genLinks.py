import json
import jsonpath_rw as jsonp
# import config.config as config
import os
import hashlib


class dataNetwork:
    # todo
    # switch to functional approach (split out the node/link matching to function,
    # map to matched nodes)
    # add 'grouping'
    # so a node becomes a group of items sharing a characteristic
    # this could be 'filename', the whole dataset ('all_papers') as now
    # or a specific mesh heading
    # searches for shared data is in or out of group
    # (i.e. all in group or all data)

    # define what to search within for linking purposes
    # i.e. 'paper' will search for associations within the current paper;
    # 'all_data' will search for associations across the whole dataset
    search_across_types = (
      'paper',
      'all_data'
    )

    def __init__(self, data_path):
        self.data_path = data_path
        self.dataset = self.getCacheData()
        self.target_path = ''
        self.search_across = self.search_across_types[0]
        self.ignore_empty_nodes = True
        self.nodes = {}
        self.links = {}

        # node_name is a list containing a json path to all the fields
        # making up the node name (used by processNodeMatch). defaults to ['$'], which will just use the
        # value of the found match. can be used to e.g. join an author's
        # first and last names from a node containing all the author data
        # as in self.node_name = ['$.family', '$.given'] will extract
        # family and given from the node value {'family':'Smith', 'given':'John'}
        self.node_name = ['$']

        # the 'additional_data' list can contain jsonpaths to bits of data we might
        # want to store with a node. Each entry is of the form
        # {'path':'[jsonpath to retrieve]', 'context':'[node|paper]', 'name':'[key name for output], 'clean_function': [cleaning_function]}
        # Specifying context as node means only this node is searched;
        # specifying context as paper searches the whole paper
        # 'name' defines the key name in each node dict
        # 'cleaning_function' allows custom cleaning function to be used
        # defaults to 'clean'
        self.additional_data_node = []

        # the below are only used if search_across is set to 'all_data'
        self.assoc_node_target_path = None
        self.link_sets = {}
        self.assoc_nodes = {}
        self.additional_data_assoc_node = []

    # get all [filenames] from cache[/filetype]
    # processes to output_type (def. json)
    # returns a dict of [filename]->[filedata]
    def getCacheData(self, output_type='json'):
        cache_files = []
        for root, dirs, files in os.walk(self.data_path):
            for name in files:
                cache_files.append(name)

        dataset = {}

        for cache_filename in cache_files:
            cache_file = open(self.data_path+'/'+cache_filename, 'r')
            if (output_type == 'json'):
                dataset[cache_filename] = json.load(cache_file)
                cache_file.close()
            else:
                raise ValueError('(genLinks.getCacheData) Error, unrecognised output_type')

        self.dataset = dataset
        return dataset

    def getPath(self, target, target_path):
        # parse the path as json
        json_path = jsonp.parse(target_path)
        # create a list of matches
        matches = [match for match in json_path.find(target)]
        # return the list of matches (DatumInContext objects, key properties: value, context (another DatumInContext object))
        return matches

    def clean(self, value):
        clean_value = value.lower().strip()
        return clean_value

    # process the found node
    # this returns the match value directly by default but may be
    # used to e.g. join multiple subfields of a match together
    def processFoundNode(self, match):
        value = ''
        if self.node_name is not None:
            for name in self.node_name:
                name_matches = self.getPath(match.value, name)
                value += ' '.join([m.value for m in name_matches])
                value += ' '
        value = self.clean(value)
        return value

    def genNetwork(self):
        nodes = {}
        assoc_nodes = {}
        complete_counter = 0
        total_items = len(list(self.dataset.items()))
        if self.search_across == 'paper':
            for paper_name, paper in list(self.dataset.items()):
                print('Generating nodes for paper {paper_name}'.format(paper_name=paper_name))
                print('Paper {n} out of {total} being processed...'.format(n=complete_counter+1, total=total_items))
                complete_counter += 1
                node_matches = self.getPath(paper, self.target_path)
                print("here0")
                if node_matches is not None and len(node_matches) > 0:
                    for match in node_matches:
                        # process the value to that set by self.node_name
                        processed_value = self.processFoundNode(match)
                        # get a cleaned value (this is what we match with)
                        clean_value = self.clean(processed_value)
                        print("here1")
                        # if we've got ignore_empty_nodes on, skip if empty string
                        if self.ignore_empty_nodes:
                            if clean_value == '':
                                continue
                        print("here")
                        node_id = hashlib.sha256(clean_value.encode('utf-8', 'ignore')).hexdigest()
                        try:
                            nodes[node_id]['count'] += 1
                        except:
                            nodes[node_id] = {
                              'value': processed_value,
                              'clean_value': clean_value,
                              'count': 1
                            }

                        # if additional_data_node is set, look for this in the relevant context
                        if len(self.additional_data_node) > 0:
                            for add_data_n in self.additional_data_node:
                                nodes[node_id][add_data_n['name']] = []
                                if add_data_n['context'] == 'node':
                                    add_data_matches = self.getPath(match.value, add_data_n['path'])
                                    for add_data_match in add_data_matches:
                                        if 'clean' in add_data_n and add_data_n['clean'] is not None:
                                            # print 'Cleaning additional_data_node {name} with function {cleaner}'.format(name=add_data_n['name'], cleaner=add_data_n['clean'].__name__)
                                            nodes[node_id][add_data_n['name']].append(add_data_n['clean'](add_data_match.value))
                                        else:
                                            nodes[node_id][add_data_n['name']].append(add_data_match.value)

                    print('Generating links for paper {paper_name}'.format(paper_name=paper_name))
                    self.getLinks(node_matches, node_id)
        elif self.search_across == 'all_data':
            # paper ($) is always root of data
            # paper has n nodes
            # paper has l link_nodes
            # gen for each paper
            # now go through each paper and match link_nodes to other papers' link_nodes
            # add to a link_sets
            # e.g. subjects relationship with other subjects when they share mesh headings
            # generate list of subjects
            # generate list of mesh headings
            # for each subject, go through each mesh heading and find in context (in this case paper). end up with count of each mesh heading per subject
            ##
            # search whole dataset to get nodes
            for paper_name, paper in list(self.dataset.items()):
                # get node(s) for this paper
                print('Generating nodes for paper {paper_name}'.format(paper_name=paper_name))
                paper_nodes = self.genNodes(paper, additional_data=self.additional_data_node)
                nodes.update(paper_nodes)

                # get association nodes for this paper
                print('Generating association nodes for paper {paper_name}'.format(paper_name=paper_name))
                paper_assoc_nodes = self.genNodes(paper, self.assoc_node_target_path, additional_data=self.additional_data_assoc_node)
                assoc_nodes.update(paper_assoc_nodes)

                # now go through nodes and assoc_nodes
                # where context contains both node and assoc_node, +1 to node-assoc_node link
                self.genNodeAssocLinks(paper_nodes, paper_assoc_nodes)
        else:
            raise ValueError('(genLinks.genNetwork) Error, search_across must be one of the values specified in search_across_types')
        self.nodes = nodes
        self.assoc_nodes = assoc_nodes
        return nodes

    def genNodes(self, target, target_path=None, additional_data=None):
        if target_path is None:
            target_path = self.target_path
        node_matches = self.getPath(target, target_path)
        nodes = {}
        if node_matches is not None and len(node_matches) > 0:
            for match in node_matches:
                processed_value = self.processFoundNode(match)
                clean_value = self.clean(processed_value)

                # if we've got ignore_empty_nodes on, skip if empty string
                if self.ignore_empty_nodes:
                    if clean_value == '':
                        continue

                node_id = hashlib.sha256(clean_value.encode('utf-8', 'ignore')).hexdigest()
                try:
                    nodes[node_id]['count'] += 1
                except:
                    nodes[node_id] = {
                      'value': processed_value,
                      'clean_value': clean_value,
                      'count': 1
                    }

                # if additional_data_node is set, look for this in the relevant context
                if len(additional_data) > 0:
                    for add_data_n in additional_data:
                        nodes[node_id][add_data_n['name']] = []
                        if add_data_n['context'] == 'node':
                            add_data_matches = self.getPath(match.value, add_data_n['path'])
                            for add_data_match in add_data_matches:
                                if 'clean' in add_data_n and add_data_n['clean'] is not None:
                                    # print 'Cleaning additional_data_node {name} with function {cleaner}'.format(name=add_data_n['name'], cleaner=add_data_n['clean'].__name__)
                                    nodes[node_id][add_data_n['name']].append(add_data_n['clean'](add_data_match.value))
                                else:
                                    nodes[node_id][add_data_n['name']].append(add_data_match.value)

        return nodes

    def genNodeAssocLinks(self, nodes, assoc_nodes):
        for node_id in nodes:
            for assoc_node_id in assoc_nodes:
                if node_id < assoc_node_id:
                    link_id = node_id+assoc_node_id
                else:
                    link_id = assoc_node_id+node_id
                link_id = hashlib.sha256(link_id).hexdigest()

                try:
                    self.links[link_id]['count'] += 1
                except:
                    self.links[link_id] = {
                      'src_id': node_id,
                      'target_id': assoc_node_id,
                      'count': 1
                    }
        return self.links

    def getLinks(self, matches, linking_node_id):
        for match in matches:
            processed_value = self.processFoundNode(match)
            clean_value = self.clean(processed_value)

            # if we've got ignore_empty_nodes on, skip if empty string
            if self.ignore_empty_nodes:
                if clean_value == '':
                    continue

            node_id = hashlib.sha256(clean_value.encode('utf-8', 'ignore')).hexdigest()
            if node_id != linking_node_id:
                if node_id < linking_node_id:
                    link_id = node_id+linking_node_id
                else:
                    link_id = linking_node_id+node_id
                link_id = hashlib.sha256(link_id).hexdigest()
                try:
                    self.links[link_id]['count'] += 1
                except:
                    self.links[link_id] = {
                      'src_id': linking_node_id,
                      'target_id': node_id,
                      'count': 1
                    }
        return self.links

    def processOutput(self):
        output = {
            'nodes': [],
            'assoc_nodes': [],
            'links': [],
        }
        for node_id in self.nodes:
            self.nodes[node_id]['id'] = node_id
            output['nodes'].append(self.nodes[node_id])
        for node_id in self.assoc_nodes:
            self.assoc_nodes[node_id]['id'] = node_id
            output['assoc_nodes'].append(self.assoc_nodes[node_id])
        for link_id in self.links:
            self.links[link_id]['id'] = link_id
            output['links'].append(self.links[link_id])
        return output
