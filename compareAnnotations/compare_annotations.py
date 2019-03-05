import sys
import os
import json
import datetime
import argparse
import re

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = '', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-g', '--genome', required = True)

args = parser.parse_args()

## MISC ########################################################################

timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

translation_locations = {'ec'      : 'KBaseOntology.OntologyTranslation.EBI_EC.ModelSEED.json',
                         'keggro'  : 'KBaseOntology.OntologyTranslation.KEGG_RXN.ModelSEED.json',
                         'keggko'  : 'KBaseOntology.OntologyTranslation.KEGG_KO.ModelSEED.json',
                         'metacyc' : 'KBaseOntology.OntologyTranslation.Metacyc_RXN.ModelSEED.json',
                         'SSO'     : 'KBaseOntology.OntologyTranslation.SSO.ModelSEED.json'}

## CLASSES #####################################################################

class Gene:
    def __init__(self, id):
        self.id = id
        self.annotations = []

    def add(self, ontology_event, type, term):
        self.annotations.append({"ontology_event" : ontology_event,
                                           "type" : type,
                                           "term" : term
                               })

class RXN:
    def __init__(self, rxn):
        self.rxn = rxn
        self.translations = []

    def add(self, ontology_event, gene, type, term):
        self.translations.append({"ontology_event" : ontology_event,
                                           "gene" : gene,
                                           "type" : type,
                                           "term" : term
                               })

## FUNCTIONS ###################################################################

def get_genome():
    genome_dict = json.loads(open(args.genome, "r").read() )
    return(genome_dict)

def get_translations(type):
    translations = json.loads(open(translation_locations[type], "r").read() )
    return(translations['translation'])

def get_translations2():
    translations = {}

    for type in translation_locations:
        ontology_translations = json.loads(open(translation_locations[type], "r").read() )
        translations[type] = {}

        for term in ontology_translations['translation']:
            for entry in ontology_translations['translation'][term]['equiv_terms']:
                rxn = entry['equiv_term']
                translations[type][term] = rxn

    return(translations)

def translate(term, translation_dict):

    translations = []

    if term in translation_dict:
        full = translation_dict[term]
        for entry in full['equiv_terms']:
            translations.append(entry['equiv_term'])
    else:
        # terms not found in the translations - collect these later
        pass

    return(translations)

def get_ontology_events(genome_dict):
    '''
    returns dict with ontology information. If description is not found, adds the method instead.
    '''

    ontology_events = {}

    for event, ontology in enumerate(genome_dict['ontology_events']):
        ontology_events[event] = {'id'          : ontology['id'],
                                  'description' : ontology.get('description', ontology['method'])
                                 }

    return ontology_events

def get_genes_and_terms2(genome_dict, genes):
    for feature in genome_dict['features']:
        gene = feature['id']
        if 'ontology_terms' in feature:
            for type in feature['ontology_terms']:
                for term in feature['ontology_terms'][type]:
                    for ontology_event in feature['ontology_terms'][type][term]:
                        if gene in genes:
                            genes[gene].add(ontology_event, type, term)
                        else:
                            genes[gene] = Gene(gene)
                            genes[gene].add(ontology_event, type, term)

    return(genes)

def get_translations(genes, rxns, genome_dict, translations, getECs = True):
    rxns['None'] = RXN('None')

    for gene in genes:
        for event in genes[gene].annotations:
            term = event['term']
            type = event['type']

            if type == 'SSO':
                term = genome_dict['ontologies_present']['SSO'][term]

            if type == 'metacyc':
                if term.startswith("META:"):
                    term = term.replace('META:', '')

            if term in translations[type]:
                rxn = translations[type][term]

                if rxn != None:
                    if rxn in rxns:
                        rxns[rxn].add(event['ontology_event'], gene, type, term)
                    else:
                        rxns[rxn] = RXN(rxn)
                        rxns[rxn].add(event['ontology_event'], gene, type, term)
                else:
                    rxns['None'].add(event['ontology_event'], gene, type, term)

            else:
                rxns['None'].add(event['ontology_event'], gene, type, term)

            if getECs: #extract ECs from SSO terms
                if type == "SSO":
                    # ECs
                    ecList = search_for_ec(term)
                    if len(ecList) > 0:
                        for ec in ecList:
                            if ec in translations['ec']:
                                rxn = translations['ec'][ec]

                                if rxn != None:
                                    if rxn in rxns:
                                        rxns[rxn].add(event['ontology_event'], gene, type, ec)
                                    else:
                                        rxns[rxn] = RXN(rxn)
                                        rxns[rxn].add(event['ontology_event'], gene, type, ec)
                                else:
                                    rxns['None'].add(event['ontology_event'], gene, type, ec)

                            else:
                                rxns['None'].add(event['ontology_event'], gene, type, ec)


    return(rxns)


def search_for_ec(line):
    ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+", line)
    return(ecList)

def summarize(genes, rxns, ontology_events):
    summary = {}

    for gene in genes:
        for event in genes[gene].annotations:
            term = event['term']
            event = event['ontology_event']

            if event not in summary:
                summary[event] = {}
            if 'gene' not in summary[event]:
                summary[event]['gene'] = []
            if 'term' not in summary[event]:
                summary[event]['term'] = []

            summary[event]['gene'].append(gene)
            summary[event]['term'].append(term)

            summary[event]['gene'] = list(set(summary[event]['gene']))
            summary[event]['term'] = list(set(summary[event]['term']))

    for rxn in rxns:
        for event in rxns[rxn].translations:
            event = event['ontology_event']
            if event not in summary:
                summary[event] = {}
            if 'rxn' not in summary[event]:
                summary[event]['rxn'] = []

            if rxn != "None":
                summary[event]['rxn'].append(rxn)

            summary[event]['rxn'] = list(set(summary[event]['rxn']))

    print("EVENT", "DESCRIPTION", "TYPE", "GENES", "TERMS", "RXNS", sep = "\t")
    for event in sorted(summary.keys()):
        description = ontology_events[event]['description']
        type = ontology_events[event]['id']
        genes_list = summary[event]['gene']
        terms_list = summary[event]['term']
        rxns_list = summary[event]['rxn']
        print(event, description, type, len(set(genes_list)), len(terms_list), len(rxns_list), sep = "\t")

    return(summary)

def make_table(rxns, ontology_events):

    table = {}

    for rxn in rxns:
        if rxn != 'None':
            table[rxn] = {}
            for event in ontology_events:
                for translation in rxns[rxn].translations:
                    if translation['ontology_event'] == event:

                        if event in table[rxn]:
                            table[rxn][event] += 1
                        else:
                            table[rxn][event] = 1
    print("RXN", end = "")
    for event in ontology_events:
        print("\t" + ontology_events[event]['description'], end = "")
    print()
    for rxn in sorted(table.keys()):
        print(rxn, end = "")
        for event in ontology_events:
            if event in table[rxn]:
                print("\t" + str(table[rxn][event]), end = "")
            else:
                print("\t0", end = "")
        print()


def convert_modelseed_to_kegg(ontology_events):
    pass

def cumulative_sum_curve2(summary, type, compare_to, ontology_events, ignore = []):

    csc = open ('csc.txt', 'w')

    working = summary.copy()
    for ignoreCheck in list(working.keys()):
        if ignoreCheck in ignore:
            print("IGNORING ONTOLOGY EVENT:", ignoreCheck)
            del working[ignoreCheck]


    cumulative_collection = [] # collection of all previously found items
    csc.write("DESCRIPTION\tADDED\tOVERLAP\tBUFFER\tTOTAL\n")

    while len(working) > 0 :

        # iterate through, removing most abundant and moving those terms to the collection list
        abundant_key = ""
        abundant_value = 0
        abundant_already_added = 0

        if compare_to == None:
            for ontology in working:

                # remove those already found in collection
                remaining = set(working[ontology][type]) - set(cumulative_collection)
                already_present =  set(cumulative_collection) & set(working[ontology][type])

                if len(remaining) > abundant_value:
                    abundant_value = len(remaining)
                    abundant_key = ontology
                    abundant_already_added = len(already_present)

        else:
            remaining = set(working[compare_to][type]) - set(cumulative_collection)
            already_present =  set(cumulative_collection) & set(working[compare_to][type])

            abundant_value = len(remaining)
            abundant_key = compare_to
            abundant_already_added = len(already_present)

            compare_to = None

        # remove most abundant
        if abundant_key == "": # if there are ontology events with 0 counts for this type
            for loser in working:
                csc.write(ontology_events[loser]['description'] + "\t" + str(0) + "\t" + str(len(working[loser][type])) + "\t" + str(len(cumulative_collection) - len(working[loser][type])) + "\t" + str(len(cumulative_collection)) + "\n")
            break

        else:
            csc.write(ontology_events[abundant_key]['description'] + "\t" +
                  str(len(working[abundant_key][type]) - abundant_already_added) + "\t" +
                  str(abundant_already_added)  + "\t" +
                  str(len(cumulative_collection) - abundant_already_added)  + "\t" +
                  str(len(cumulative_collection) + len(working[abundant_key][type]) - abundant_already_added) + "\n")

            cumulative_collection = set(cumulative_collection) | set(working[abundant_key][type])
            working.pop(abundant_key)

    csc.close()

    csc_filename = args.genome + "_" + type + ".png"
    csc_title = args.genome + "_" + type

    os.system('Rscript csc.R csc.txt ' + csc_filename + " " + csc_title)
    print("\n*** Cumulative sum plot data written to csc.txt and plot written to " + csc_filename + "\n")

def calculate_overlaps(summary, type):
    items = {}
    for ontology in summary:
        for item in summary[ontology][type]:
            if item in items:
                items[item].append(ontology)
            else:
                items[item] = [ontology]

    combinations = {}

    for item in items:
        #print(item, sorted(items[item]), sep = "\t")
        if str(sorted(items[item])) in combinations:
            combinations[str(sorted(items[item]))] += 1
        else:
            combinations[str(sorted(items[item]))] = 1

    for combo in combinations:
        print(combo, combinations[combo], sep = "\t")

################################################################################

def main():

    # Prepare Data
    genome_dict = get_genome()
    ontology_events = get_ontology_events(genome_dict)
    translations = get_translations2()
    genes = {} # holds instances of the Gene class
    rxns = {}  # holds instances of the RXN class
    genes = get_genes_and_terms2(genome_dict, genes)
    rxns = get_translations(genes, rxns, genome_dict, translations, getECs = True)

    summary = summarize(genes, rxns, ontology_events)
    cumulative_sum_curve2(summary, 'gene', 0, ontology_events, ignore = [])
    #calculate_overlaps(summary, 'rxn')
    #make_table(rxns, ontology_events)
    #print(ontology_events)





## RUN #########################################################################

main()
