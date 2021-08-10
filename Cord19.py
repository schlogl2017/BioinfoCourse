#!/usr/bin/env python
# coding: utf-8

# In[12]:


import os
import glob
import re
import csv
import functools
import operator
# import itertools
import json
# import nltk
# from nltk.tokenize import sent_tokenize
import pandas as pd
from collections import defaultdict
from tqdm.notebook import tqdm
from pandarallel import pandarallel
# import matplotlib.pyplot as plt

# plt.style.use('ggplot')


# In[2]:


pandarallel.initialize(use_memory_fs=False,nb_workers=2)


# In[2]:


get_ipython().system('ls Data/archive/')


# In[3]:


data = pd.read_csv('Data/archive/metadata.csv',
                  usecols=['cord_uid','source_x','title','license','publish_time',
                           'abstract','authors','journal','url'])
data.head()


# In[4]:


data.tail()


# In[3]:


dirs = ["Data/archive/document_parses/pdf_json", "Data/archive/document_parses/pmc_json"]


# In[4]:


all_jsons = (glob.glob(di+'/*.json') for di in dirs)


# In[5]:


# all_files = itertools.chain.from_iterable(all_jsons)

all_files = functools.reduce(operator.iconcat, all_jsons, [])


# In[57]:


all_files[1]


# In[6]:


len(all_files)


# In[11]:


pd.DataFrame.from_dict({"doc_id": [None], 
                        "source": [None], 
                        "title": [None],
                        "abstract": [None], 
                        "text_body": [None]})
        


# In[6]:


def get_authors(json_file_metadata_authors):
    """
    Inputs:
        json_file_metadata_aoutors = list
    """
    names = []
    for a in json_file_metadata_authors:
        first =  a['first']
        last = a['last']
        middle = ''.join(a['middle'])
        if middle != '':
            names.append((first, middle, last))
        else:
            names.append((first, last))
    return [' '.join(name) for name in names if '' not in name]


# In[50]:


get_authors(j['metadata']['authors'])


# In[22]:


def format_name(author):
    middle_name = " ".join(author['middle'])
    
    if author['middle']:
        return " ".join([author['first'], middle_name, author['last']])
    else:
        return " ".join([author['first'], author['last']])


# In[36]:


j['metadata']['authors'][3]


# In[37]:


[(di['section'], di['text']) for di in j['body_text']]


# In[38]:


{di['section']: "" for di in j['body_text']}


# In[7]:


def format_body(body_text):
    texts = [(di['section'], di['text']) for di in body_text]
    texts_di = {di['section']: "" for di in body_text}
    
    for section, text in texts:
        texts_di[section] += text

    body = ""

    for section, text in texts_di.items():
        #body += section
        #body += "\n\n"
        body += text
        #body += "\n\n"
    
    return body


# In[40]:


format_body(j['body_text'])


# In[46]:


text == format_body(j['body_text'])


# In[65]:


j['metadata']["title"]


# In[285]:


# def parse_body_text(body_text):
#     body =""
#     for item in body_text:
#         body += item["section"]
#         body += "\n\n"
#         body += item["text"]
#         body += "\n\n"
#     body = clean(body)
#     return body


# In[8]:


def subs_months(month):
    month_map = {
        'January':'1',
        'February':'2',
        'March':'3',
        'April':'4',
        'May':'5',
        'June':'6',
        'July':'7',
        'August':'8',
        'September':'9',
        'October':'10',
        'November':'11',
        'December':'12'
    }
    m_splt = month.split('/')[1]
    new = month.replace(m_splt, month_map[m_splt])
    return new


#     1 - \[.*?\] = [{'start': 215, 'end': 216, 'text': '8', 'ref_id': None}], [],[17] [18] [19]...
#     2 - \(.*?\) = (COVID- 19), (in China and Russia, specifically), (RBD), (sfGFP), (ACE2), (Fig. 1a), ...
#     3 - "\s+ = \s pontons/matches any whitespace character (equivalent to [\r\n\t\f\v  ]) + matches the previous token between one and unlimited times, as many times as possible, giving back as needed (greedy)
#     4 - \w*\d\w* = \w matches any word character (equivalent to [a-zA-Z0-9_]) * matches the previous token between zero and unlimited times, as many times as possible, giving back as needed (greedy) \d matches a digit (equivalent to [0-9])
#     5 - \w+…|… = 1st Alternative \w+…\w matches any word character (equivalent to [a-zA-Z0-9_]) + matches the previous token between one and unlimited times, as many times as possible, giving back as needed (greedy) … matches the character … with index 823010 (202616 or 200468) literally (case sensitive) 2nd Alternative …… matches the character … with index 823010 (202616 or 200468) literally (case sensitive)
# 

# In[9]:


def clean(text):
    punctuation = '!"#$%&\'()*+,/:;<=>?@[\\]^_`{|}~'
    text = text.lower()
    #clean brackts and its content
    text = re.sub(r'\[.*?\]', '', text)
    # clean all inside ( )
    text = re.sub(r'\(.*?\)', '', text)
    # matchs spaces and the preceding character
    text = re.sub(r"\s+", " ", text)
    text = re.sub(r"\w+':\s", '', text)
    text = text.replace('s ars-cov-2', 'sars-cov-2')
    text = re.sub(r"\w+…|…", "", text)  # Remove ellipsis (and last word)
    text = re.sub(f"[{re.escape(punctuation)}]", " ", text)
    return text


# In[53]:


key_list = ["paper_id", "paper_title", "authors", "abstract", "body_text"]
json_data = defaultdict(list, [(k, [None]) for k in key_list])


# In[73]:


j['metadata']['title']


# In[93]:


for entry in j['body_text']:
    print(entry['text'])


# In[10]:


def get_json_data(file):
    key_list = ["paper_id", "paper_title", "authors", "abstract", "body_text"]
    json_data = defaultdict(str, [(k, None) for k in key_list])
    with open(file, 'r') as fh:
        js = json.load(fh)
        json_data["paper_id"] = js['paper_id']
        try:
            for entry in js['abstract']:
                json_data["abstract"] = entry['text']
        except:
            pass
        try:
            text =""
            for entry in js['body_text']:
                text += entry['text']
            json_data["body_text"] = text
        except:
            pass
        try:
            for entry in js['metadata']:
                json_data["paper_title"] = js['metadata']['title']
                json_data["authors"] = get_authors(js['metadata']['authors'])
        except:
            pass
    return json_data


# In[121]:


get_json_data('Data/archive/document_parses/pdf_json/b8e69fa00cb9e90e48340337cf69db7dc76db884.json')


# In[113]:


# def get_data(file):
#     paper_id = ""
#     paper_title = ""
#     authors = []
#     abstract = ""
#     body_text = ""
#     date = ""
#     doi = ""
#     with open(file, 'r') as fh:
#         js = json.load(fh)
#         paper_id += js['paper_id']
#         paper_title += js['metadata']['title']
#         authors = get_authors(js['metadata']['authors'])
#         publisher = " "
#         try:
#             for entry in js['abstract']:
#                 abstract += entry['text']
#         except:
#             abstract += ""
#         try:
#             for entry in js['body_text']:
#                 body_text += entry['text']
#         except:
#             body_text += ""
#         try:
#             d = '/'.join(js['body_text'][-1]['text'].split(';')[1].split(' ')[2:])
#             date += subs_months(d)
#             publisher += js['body_text'][24]['text'].split(' ')[11]
#             doi += ''.join(js['body_text'][-2]['text'].split(' ')[-2:])
#         except:
#             date += ""
#             doi += ""
#             publisher += ""
#     return paper_id, paper_title, authors, abstract, body_text


# In[13]:


j = json.load(open('Data/archive/document_parses/pdf_json/b8e69fa00cb9e90e48340337cf69db7dc76db884.json'))


# In[44]:


Id, tit, aut, abst, text = get_data('Data/archive/document_parses/pdf_json/b8e69fa00cb9e90e48340337cf69db7dc76db884.json')


# In[152]:


d


# In[86]:


r = get_data('Data/archive/document_parses/pdf_json/fc0b4e76a1dc550b3c147ce8cc92516b6364a8fa.json')


# In[14]:


for key in j.keys(): print(key)


# In[111]:


len(j['ref_entries'])


# In[112]:


for i, d in enumerate(j['ref_entries']):
    print(i, d)
    


# In[149]:


j['body_text'][24]['text'].split(' ')[11]


# In[ ]:


publisher = js['body_text'][24]['text'].split(' ')[11]


# In[120]:


doi = ''.join(js['body_text'][-2]['text'].split(' ')[-2:])


# In[147]:


j['body_text'][-1]['text']


# In[130]:


date = '/'.join(j['body_text'][-1]['text'].split(';')[1].split(' ')[2:])


# In[131]:


date


# In[133]:


subs_months(date)


# In[69]:


clean(body_text)


# In[156]:


text_data = (get_data(file) for file in ['Data/archive/document_parses/pdf_json/b8e69fa00cb9e90e48340337cf69db7dc76db884.json',
                                         'Data/archive/document_parses/pdf_json/fc0b4e76a1dc550b3c147ce8cc92516b6364a8fa.json'])


# In[157]:


pd.DataFrame(text_data, columns=[
    "paper_id", "paper_title", "authors", "abstract", "body_text", "publisher", "date"
])


# In[241]:


# t = ' '.join([str(j_file['body_text'][i]) for i in range(len(text))])


# In[234]:


# def json_authors_clean(json_file_metadata_authors):
#     """
#     Inputs:
#         json_file_metadata_authors = list
#     """
#     names = []
#     for a in json_file_metadata_authors:
#         first =  a['first']
#         last = a['last']
#         middle = ''.join(a['middle'])
#         if middle != '':
#             names.append((first, middle, last))
#         else:
#             names.append((first, last))
#     return [' '.join(name) for name in names if '' not in name]            


# In[114]:


text_data = (get_data(file) for file in all_files)


# In[19]:


data = (get_json_data(file) for file in all_files)


# In[ ]:


csv_columns = ['paper_id', 'paper_title', 'authors', 'abstract', 'body_text']

csv_file = 'text_file_from_json.csv'
try:
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for doc in data:
            writer.writerow(doc)
except IOError:
    print("I/O error")


# In[27]:


#dict from csv
# mydict = {y[0]: y[1] for y in [x.split(",") for x in open('file.csv').read().split('\n') if x]}

#or
# input_file = csv.DictReader(open("coors.csv"))


# In[5]:


csv_file = "text_file_from_json.csv"


# In[170]:


get_ipython().system('cat text_file_from_json.csv | head -n100')


# In[30]:


get_ipython().run_line_magic('pinfo', 'pd.read_csv')


# In[116]:


names=["paper_id", "paper_title", "authors", "abstract", "body_text"]


# In[ ]:


df = pd.DataFrame.from_records(data, columns=names)


# In[ ]:


def search_dataframe(df,search_words):
    search_words=stem_words(search_words)
    df1=df[functools.reduce(lambda a, b: a&b, (df['abstract'].str.contains(s) for s in search_words))]
    return df1


# In[7]:


df.head()


# In[8]:


df['paper_title'].isnull().sum()


# In[9]:


df['paper_id'].isnull().sum()


# In[10]:


df['authors'].isnull().sum()


# In[11]:


df['abstract'].isnull().sum()


# In[12]:


df['body_text'].isnull().sum()


# In[13]:


df['publisher'].isnull().sum()


# In[14]:


df['date'].isnull().sum()


# In[16]:


df['date'].unique()


# In[ ]:




