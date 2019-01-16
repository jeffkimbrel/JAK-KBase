import sys
import os
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = '', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-k', '--kegg',
    required = True)

args = parser.parse_args()

## MAKE MAP ###

map = {}
map_call = 'curl -g -s -S http://rest.kegg.jp/link/reaction/orthology/'
map_raw = os.popen(map_call).read()
for line in map_raw.split("\n"):
    split = line.split("\t")

    if len(split) > 1:
        ko = split[0].replace("ko:", "")
        rxn = split[1].replace("rn:", "")

        if ko in map:
            map[ko].append(rxn)
        else:
            map[ko] = [rxn]

## MAP STUFF ##

kegg_raw = [line.strip() for line in open(args.kegg)]
for line in kegg_raw:
    split = line.split("\t")
    if len(split) > 1:

        if split[1] in map:
            reactions = map[split[1]]
            for rn in reactions:
                print(split[0], rn, sep = "\t")
