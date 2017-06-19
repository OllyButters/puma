from pyzotero import zotero
import json
import config.config as config
import re

# extends the pyzotero Zotero class
# -
# adds a setter to retrieve the collection key (if present) when setting collection as a human-readable name
# -
# adds method to retrieve all the keys of items in the collection or repo
# adds method to retrieve all the items in the collection or repo
class zotPaper (zotero.Zotero):

  def __init__(self, **kwargs):
    self.__collection = None
    self.collection_key = None
    self.papers = []
    super(zotPaper, self).__init__(config.zotero_id, config.zotero_type, config.zotero_api_key, **kwargs)

  @property
  def collection(self):
    return self.__collection

  @collection.setter
  def collection(self, collection):
    self.__collection = collection
    if collection is not None:
      self.collection_data = self.collections()
      for collection in self.collection_data:
        if collection['data']['name'] == self.__collection:
          self.collection_key = collection['data']['key']
          break
      else:
        self.collection_key = None
    else:
      self.collection_key = None

  # get a list of item keys for this instance
  def getPapersKeys(self, *args, **kwargs):
    # get all keys
    if self.collection_key is not None:
      zot_keys = self.collection_items(format='keys', limit=999999)
    else:
      zot_keys = self.items(format='keys', limit=999999)
    self.papers_keys = zot_keys.split('\n')
    return self.papers_keys

  # populate self.papers list with all papers in self.collection
  def getPapersList(self, key_list = None, *args, **kwargs):
    self.papers = []

    if key_list is None:
      # get all keys
      key_list = self.getPapersKeys()

    # get items one by one
    # if we just use the items() method, the not all data is retrived
    # we therefore request individually by key
    for key in key_list:
      if key != '':
        print 'Getting zotero item '+key
        item = super(zotPaper, self).item(key, *args, **kwargs)
        self.papers.append(item)
    return self.papers


  def extraToFields(self, extra):
    entries = extra.split("\n");
    data = {}
    for entry in entries:
      datum = re.search(r'^(.+?):(.+?)$', entry)
      if datum is not None:
        try:
          data[datum.group(1)] = datum.group(2)
        except IndexError:
          #if either value is missing just add to the 'extra' value as text
          data['extra'] += entry+"\n"
      else:
        data['extra'] += entry+"\n"
        # don't raise an error, this will cause any other unformatted extra
        # to vanish
    return data


  def fieldsToExtra(self, paper, fields):
    extra = paper['extra']
    for field in fields:
      value = ''
      if field in paper:
        value = paper[field]
      extra += '%s:%s\n' % (field, paper[field])
    return extra


  def uploadExtra(self, paper):
    extra = paper['extra']
    extra_dict = self.extraToFields(extra)
    key = paper['key']
    try:
      # get the current data
      zot_paper = self.item(key)
      zot_extra_dict = {}
      if 'extra' in zot_paper['data']:
        zot_extra_dict = self.extraToFields(zot_paper['data']['extra'])
    except:
      # do logging here
      print "Error getting paper %s from Zotero" % (key)

    # compare values and update zot_extra_dict only if key not present
    for k,v in extra_dict.iteritems():
      if k not in zot_extra_dict.keys():
        zot_extra_dict[k] = v

    # now flatten zot_extra_dict ready for upload
    zot_extra = ''
    for k,v in zot_extra_dict.iteritems():
      zot_extra += '%s="%s"\n' % (k,v)

    zot_paper['data']['extra'] = zot_extra

    self.update_item(zot_paper)
    return zot_paper

  # map keys in paper to zotero dict
  # defines mappings as dicts by src_type, [src_type key]: [zotero key]
  # used for quick conversions of lists of papers (as dicts) to zotero format ready for upload
  # for complex conversions use the Merge class in papersMerge
  def mapFields(self, paper, item_type = 'journalArticle', src_type = 'doi'):
    fields = self.item_type_fields(item_type)
    creator_types = self.item_creator_types(item_type)
    creator_fields = self.creator_fields()
    all_fields = self.item_fields()

    creator_mappings = {
      'doi': {
        'fieldname': 'author',
        'fields': {
          'given': 'firstName',
          'family': 'lastName',
          'name': 'name',
        }
      }
    }

    field_mappings = {
      'doi': {
        'title': 'title',
        'DOI': 'DOI',
        'ISSN': 'issn',
        'container-title': 'publicationTitle',
        'volume': 'volume',
        'pages': 'page',
        'URL': 'url',
        'issue': 'issue',
        'notes': 'extra',
      }
    }

    if src_type in field_mappings:
      field_map = field_mappings[src_type]
    else:
      raise ValueError('(zotPaper::mapFields) src_type not recognised in field_mappings')

    if src_type in creator_mappings:
      creator_field = creator_mappings[src_type]['fieldname']
      creator_map = creator_mappings[src_type]['fields']
    else:
      raise ValueError('(zotPaper::mapFields) src_type not recognised in creator_mappings')

    zot_paper = {}
    zot_paper['itemType'] = item_type

    for field in paper:
      if field in field_map.keys():
        #hacky fix for issn
        if field == 'issn' and isinstance(paper[field], list):
          zot_paper[field_map[field]] = ','.join(paper[field])
        zot_paper[field_map[field]] = paper[field]
      elif field == creator_field:
        authors = paper[field]
        creators = []
        for author in authors:
          creator = {}
          creator['creatorType'] = 'author'
          for a_field in creator_map:
            if a_field in author.keys():
              creator[creator_map[a_field]] = author[a_field]
          creators.append(creator)
        zot_paper['creators'] = creators

    return zot_paper
