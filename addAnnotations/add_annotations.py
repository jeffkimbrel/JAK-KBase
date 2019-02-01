import sys
import os
import json
import datetime
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = '', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-n', '--namespace', required = True)
parser.add_argument('-o', '--out', required = True)
parser.add_argument('-g', '--genome', required = True)
parser.add_argument('-a', '--annotations', required = True)

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

            else:
                name = ""

            # adds valid and not valid annotations for recordkeeping
            ontologyCheck = {"id"    : id,
                             "name"  : name,
                             "valid" : valid
                            }

            self.ontologyChecked.append(ontologyCheck)

    def hasValidAnnotations(self):
        valid = False
        if len(self.ontologyChecked) > 0:
            for annotation in self.ontologyChecked:
                if annotation['valid'] == 1:
                    valid = True
        return(valid)


## FUNCTIONS ###################################################################

def genome_to_dict():
    genome_dict = json.loads(open(args.genome, "r").read() )
    return(genome_dict)

def ontology_dictionary_to_dict():
    ontology_dict = {}

    # KO
    if args.namespace == "KO":
        ontology_dict_raw = json.loads(open("KEGG_KO_ontologyDictionary.json", "r").read() )
    elif args.namespace == "RXN":
        ontology_dict_raw = json.loads(open("KEGG_RXN_ontologyDictionary.json", "r").read() )

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

def add_ontology_event(genome_dict):
    genome_dict['ontology_events'].append(
        {
            "id"             : args.namespace,
            "method"         : "TEST",
            "method_version" : "TEST",
            "ontology_ref"   : "TEST",
            "timestamp"      : timestamp
        }
    )

    return(genome_dict)

def summarize(genes):

    validGeneCount = 0
    invalidGenes = []
    validOntologyTermCount = 0
    invalidOntologyTerms = []

    for gene in genes:
        if genes[gene].valid == 1:
            validGeneCount += 1
        elif genes[gene].valid == 0:
            invalidGenes.append(genes[gene].id)

        for annotation in genes[gene].ontologyChecked:
            if annotation['valid'] == 1:
                validOntologyTermCount += 1
            elif annotation['valid'] == 0:
                invalidOntologyTerms.append(annotation['id'])

    print(validGeneCount, len(invalidGenes), str(invalidGenes), sep = "\t")
    print(validOntologyTermCount, len(invalidOntologyTerms), str(invalidOntologyTerms), sep = "\t")

## MAIN ########################################################################

def main():

    # read genome and ontology dictionary into to their own dictionaries
    genome_dict = genome_to_dict()
    ontology_dict = ontology_dictionary_to_dict()

    genome_dict = add_ontology_event(genome_dict)
    current_ontology_event = len(genome_dict['ontology_events']) - 1

    # read annotations and create non-validated gene classes
    annotations_to_genes()

    for gene in genes:
        genes[gene].validateGeneID(genome_dict)
        genes[gene].validateAnnotationID(ontology_dict)

    # add to genome dict
    for feature in genome_dict['features']:

        geneID = feature['id']

        if geneID in genes:
            if genes[geneID].hasValidAnnotations() == True:

                # create some things if they don't exist
                if 'ontology_terms' not in feature:
                    feature['ontology_terms'] = {}

                if args.namespace not in feature['ontology_terms']:
                    feature['ontology_terms'][args.namespace] = {}

                for annotation in genes[geneID].ontologyChecked:
                    if annotation['valid'] == 1:

                        # add to ontologies present
                        if args.namespace not in genome_dict['ontologies_present']:
                            genome_dict['ontologies_present'][args.namespace] = {}

                        if annotation['id'] not in genome_dict['ontologies_present'][args.namespace]:
                            genome_dict['ontologies_present'][args.namespace][annotation['id']] = annotation['name']

                        if annotation['id'] not in feature['ontology_terms'][args.namespace]:
                            feature['ontology_terms'][args.namespace][annotation['id']] = [current_ontology_event]
                        else:
                            feature['ontology_terms'][args.namespace][annotation['id']].append(current_ontology_event)

    with open(args.out, 'w') as outfile:
        json.dump(genome_dict, outfile, indent = 2)

    # sort ontologies present
    for ontology in genome_dict['ontologies_present']:
        genome_dict['ontologies_present'][ontology] = sorted(genome_dict['ontologies_present'][ontology].keys())

    for pos, event in enumerate(genome_dict['ontology_events']):
        print(pos, event)

    summarize(genes)

main()
