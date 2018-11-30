import sys
import os
import json
import datetime

timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

kegg_dict = {}

# get kegg release

kegg_api_release = 'curl -g -s -S http://rest.kegg.jp/info/kegg'
kegg_release_raw = os.popen(kegg_api_release).read()
kegg_release = kegg_release_raw.split("\n")[1].split()[2]

kegg_dict = {'data_version'  : kegg_release,
            'date'           : timestamp,
            'format_version' : 'N/A',
            'ontology'       : 'kegg_orthology'
            }

# get KEGG data

'''
for testing, use the dowloaded data... otherwise, on comment these two lines and
remove the call to the file
'''

#kegg_api_call = 'curl -g -s -S http://rest.kegg.jp/list/orthology'
#kegg_raw = os.popen(kegg_api_call).read()

with open('kegg_orthologs.txt', 'r') as myfile:
  kegg_raw = myfile.read()

'''

'''

kegg_dict['term_hash'] = {}

for line in kegg_raw.split("\n"):

    if len(line) > 0:

        id, name_raw = line.split("\t")
        id = id.replace("ko:", "")

        nameSplit = name_raw.split(";")

        name = ""
        synonym = ""

        if len(nameSplit) > 1:
            name = nameSplit[1].strip()
            synonym = nameSplit[0].strip()
        else:
            print(name)
            #name = nameSplit

        kegg_dict['term_hash'][id] = {'id' : id, 'name' : name, 'synonyms' : synonym.split(", ")}

with open('test.json', 'w') as outfile:
    json.dump(kegg_dict, outfile, indent = 2)
