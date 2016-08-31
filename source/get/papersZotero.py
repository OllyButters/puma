from pyzotero import zotero
import json
import config.config as config
import re

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


  #populate self.papers list with all papers in self.collection
  def getPapersList(self, *args, **kwargs):
    self.papers = []
    #get num items in collection
    if self.collection_key is not None:
      num_items = self.num_collectionitems(self.collection_key)
    else:
      num_items = self.num_items()
    #TEMP SOLUTION todo use collection_items(format="keys") to get list of keys from collection and then call item() on each new
    #as collection_items returns a limit of 100 items we need make several requests to get all the items in the collection
    for i in range(0, num_items, 100):
      if self.collection_key is not None:
        subset = super(zotPaper, self).collection_items(self.collection_key, start=i, *args, **kwargs)
      else:
        subset = super(zotPaper, self).items(start=i, *args, **kwargs)

      self.papers = self.papers + subset
    return self.papers

  def notesToFields(self, notes):
    entries = notes.split("\n");
    data = {}
    for entry in entries:
      datum = re.search(r'^(.+?)="(.+?)"$', entry)
      if datum is not None:
        data[datum.group(1)] = datum.group(2)
      else:
        raise ValueError('datum is not formatted correctly')
    return data

  def fieldsToNotes(self, paper, fields):
    notes = ''
    for field in fields:
      value = ''
      if field in paper:
        value = paper[field]
      notes += '%s="%s"\n' % (field, paper[field])
    return notes
      
  def uploadNotes(self, paper):
    notes = paper['notes']
    notes_dict = self.notesToFields(notes)
    key = paper['key']
    try:
      #get the current data
      zot_paper = self.item(key)
      zot_notes_dict = {}
      if 'notes' in zot_paper['data']:
        zot_notes_dict = self.notesToFields(zot_paper['data']['notes'])
    except:
      #do logging here
      print "Error getting paper %s from Zotero" % (key)
    
    #compare values and update zot_notes_dict only if key not present
    for k,v in notes_dict.iteritems():
      if k not in zot_notes_dict.keys():
        zot_notes_dict[k] = v

    #now flatten zot_notes_dict ready for upload
    zot_notes = ''
    for k,v in zot_notes_dict.iteritems():
      zot_notes += '%s="%s"\n' % (k,v)

    zot_paper['data']['notes'] = zot_notes

    self.update_item(zot_paper)
    return zot_paper
      
  def mapFields(self, paper, item_type = 'journalArticle', creator_type = 'author', src_type = 'doi'):
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
      if field in field_map:
        zot_paper[field_map[field]] = paper[field]
      elif field == creator_field:
        authors = paper[field]
        creators = []
        for author in authors:
          creator = {}
          creator['creatorType'] = 'author'
          for a_field in creator_map:
            if a_field in author:
              creator[creator_map[a_field]] = author[a_field]
          creators.append(creator)
        zot_paper['creators'] = creators

    return zot_paper

  


        
        
  
  
