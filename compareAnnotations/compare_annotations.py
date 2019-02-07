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

translation_locations = {'ec'     : 'KBaseOntology.OntologyTranslation.EBI_EC.ModelSEED.json',
                         'keggro' : 'KBaseOntology.OntologyTranslation.KEGG_RXN.ModelSEED.json',
                         'metacyc' : 'KBaseOntology.OntologyTranslation.Metacyc_RXN.ModelSEED.json'}

## FUNCTIONS ###################################################################

def get_genome():
    genome_dict = json.loads(open(args.genome, "r").read() )
    return(genome_dict)

def get_ec_to_modelseed():
    ec_to_modelseed = json.loads(open(translation_locations['ec'], "r").read() )
    return(ec_to_modelseed['translation'])

def get_keggro_to_modelseed():
    keggro_to_modelseed = json.loads(open(translation_locations['keggro'], "r").read() )
    return(keggro_to_modelseed['translation'])

def get_metacyc_to_modelseed():
    metacyc_to_modelseed = json.loads(open(translation_locations['metacyc'], "r").read() )
    return(metacyc_to_modelseed['translation'])

def translate(term, translation_dict):

    translations = []
    if term in translation_dict:
        full = translation_dict[term]
        for entry in full['equiv_terms']:
            translations.append(entry['equiv_term'])

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

def get_genes_and_terms(genome_dict, ontology_events):

    for ontology in ontology_events:
        ontology_events[ontology]['genes'] = []
        ontology_events[ontology]['terms'] = []

    for feature in genome_dict['features']:
        if 'ontology_terms' in feature:
            for term in feature['ontology_terms']:
                for entry in feature['ontology_terms'][term]:
                    for ont in feature['ontology_terms'][term][entry]:
                        if feature['id'] not in ontology_events[ont]['genes']:
                            ontology_events[ont]['genes'].append(feature['id'])
                        if entry not in ontology_events[ont]['terms']:
                            ontology_events[ont]['terms'].append(entry)

    return ontology_events

def search_for_ec(line):
    ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]", line)
    return(ecList)

def convert_terms_to_modelseed(genome_dict, ontology_events):
    for ontology in ontology_events:
        ontology_events[ontology]['modelseed'] = []

        if ontology_events[ontology]['id'] == 'ec':
            ec_to_modelseed = get_ec_to_modelseed()
            for term in ontology_events[ontology]['terms']:
                translations = translate(term, ec_to_modelseed)
                ontology_events[ontology]['modelseed'] += translations

        elif ontology_events[ontology]['id'] == 'keggko':
            pass

        elif ontology_events[ontology]['id'] == 'metacyc':
            metacyc_to_modelseed = get_metacyc_to_modelseed()

            for term in ontology_events[ontology]['terms']:
                if term.startswith("META:"):
                    term = term.replace('META:', '')

                translations = translate(term, metacyc_to_modelseed)
                ontology_events[ontology]['modelseed'] += translations

        elif ontology_events[ontology]['id'] == 'keggro':
            keggro_to_modelseed = get_keggro_to_modelseed()

            for term in ontology_events[ontology]['terms']:
                translations = translate(term, keggro_to_modelseed)
                ontology_events[ontology]['modelseed'] += translations

        elif ontology_events[ontology]['id'] == 'SSO':
            ec_to_modelseed = get_ec_to_modelseed()
            sso_dict = genome_dict['ontologies_present']['SSO']
            for term in ontology_events[ontology]['terms']:
                if sso_dict[term] != 'Unknown':
                    ecList = search_for_ec(sso_dict[term])
                    if len(ecList) > 0:
                        for term in ecList:
                            translations = translate(term, ec_to_modelseed)
                            ontology_events[ontology]['modelseed'] += translations

        ontology_events[ontology]['modelseed'] = sorted(list(set(filter(None, ontology_events[ontology]['modelseed']))))

    return(ontology_events)

def summarize_ontology_events(ontology_events):
    print("DESCRIPTION", "GENES", "TERMS", "RXNs", sep = "\t")
    for ontology in ontology_events:
        print(ontology_events[ontology]['description'], len(ontology_events[ontology]['genes']), len(ontology_events[ontology]['terms']), len(ontology_events[ontology]['modelseed']), sep = "\t")

def cumulative_sum_curve(ontology_events, type, compare_to):

    csc = open ('csc.txt', 'w')

    working = ontology_events.copy()
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
                csc.write(working[loser]['description'] + "\t" + str(0) + "\t" + str(len(working[loser][type])) + "\t" + str(len(cumulative_collection) - len(working[loser][type])) + "\t" + str(len(cumulative_collection)) + "\n")
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

def calculate_overlaps(ontology_events, type):
    print(ontology_events)

def main():
    genome_dict = get_genome()
    ontology_events = get_ontology_events(genome_dict)
    ontology_events = get_genes_and_terms(genome_dict, ontology_events)
    ontology_events = convert_terms_to_modelseed(genome_dict, ontology_events)

    # make comparisons
    summarize_ontology_events(ontology_events)

    cumulative_sum_curve(ontology_events, 'terms', 0)
    #calculate_overlaps(ontology_events, 'modelseed')

## RUN

main()
