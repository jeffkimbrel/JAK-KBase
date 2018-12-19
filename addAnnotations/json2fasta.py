import sys
import json

jsonFile = json.loads(open(sys.argv[1], "r").read() )
for feature in jsonFile['features']:
    print(">" + feature['id'] + "\n" + feature['protein_translation'])
