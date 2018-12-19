import sys
import os
import json
import datetime
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = '', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-n', '--namespace',
    default = "KO",
    required = False)

parser.add_argument('-g', '--genome',
    required = True)

parser.add_argument('-a', '--annotations',
    required = True)

args = parser.parse_args()


## MISC ########################################################################

timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


## CLASSES #####################################################################

genes = {}
class Gene:
    def __init__(self, id):
        self.id = id
        self.valid = 0
        self.annotations = []
        self.ontologyChecked = []

    def addAnnotation(self, annotation):
        self.annotations.append(annotation)
        self.annotations = list(set(self.annotations))

    def validateGeneID(self, genome_dict):
        for feature in genome_dict['features']:
            if feature['id'] == self.id:
                self.valid = 1

    def validateAnnotationID(self, ontology_dict):
        for id in self.annotations:
            valid = 0
            if id in ontology_dict:
                valid = 1
                name = ontology_dict[id]
                print(id, valid, name, sep = "\t")
            else:
                name = ""

            ontologyCheck = {"id"    : id,
                             "name"  : name,
                             "valid" : valid
                            }

            self.ontologyChecked.append(ontologyCheck)


## FUNCTIONS ###################################################################

def genome_to_dict():
    genome_dict = json.loads(open(args.genome, "r").read() )
    return(genome_dict)

def ontology_dictionary_to_dict():
    ontology_dict = {}

    # KO
    if args.namespace == "KO":
        ontology_dict_raw = json.loads(open("KEGG_KO_ontologyDictionary.json", "r").read() )

    # !---> add other ontology namespaces here <---! #

    # process into dict
    for entry in ontology_dict_raw['term_hash']:

        id = ontology_dict_raw['term_hash'][entry]['id']
        name = ontology_dict_raw['term_hash'][entry]['name']
        ontology_dict[id] = name

    return(ontology_dict)

def annotations_to_genes():
    annotations_raw = [line.strip() for line in open(args.annotations)]
    for line in annotations_raw:
        if not line.startswith('#'): # ignore comment lines
            elements = line.split("\t") # can add commas here as well for .csv

            if elements[0] != "": # ignore blank lines
                geneID = elements[0]
                annotation = ""

                # make Gene class for geneID
                if geneID not in genes:
                    genes[geneID] = Gene(geneID)

                # and add the (not yet validated) annotations, columns above 2 are ignored
                if len(elements) > 1:
                    annotation = elements[1]
                    genes[geneID].addAnnotation(annotation)


## MAIN ########################################################################

def main():

    # read genome and ontology dictionary into to their own dictionaries
    genome_dict = genome_to_dict()
    ontology_dict = ontology_dictionary_to_dict()

    # read annotations and create non-validated gene classes
    annotations_to_genes()

    for gene in genes:
        print("\n---", gene, "---")
        print(genes[gene].__dict__)

        genes[gene].validateGeneID(genome_dict)
        genes[gene].validateAnnotationID(ontology_dict)

        print(genes[gene].__dict__)


main()
