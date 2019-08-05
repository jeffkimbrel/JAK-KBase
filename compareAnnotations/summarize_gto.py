import json


def summarize(gto, translations):
    summary = {"genes": {},
               "terms": {},
               "rxns": {},
               "ontology_events": {},
               "orphan_terms": {}
               }

    # add ontology events
    for count, oe in enumerate(gto['ontology_events']):
        summary['ontology_events'][count] = oe

    # add gene id to summary
    for feature in gto['features']:
        gene_id = feature['id']
        summary["genes"][gene_id] = {"terms": {},
                                     "rxns": {}
                                     }

        # get ontology term
        if "ontology_terms" in feature:
            for type in feature['ontology_terms']:
                term_dict = feature['ontology_terms'][type]

                for term in term_dict:
                    for oe in term_dict[term]:

                        rxn = "none"

                        # get rxn
                        ontology_type = summary['ontology_events'][oe]['id']

                        # fix metacyc terms
                        if ontology_type == 'metacyc':
                            if term.startswith("META:"):
                                term = term.replace('META:', '')

                        # fix SSO terms
                        if ontology_type == 'SSO':

                            if term in gto['ontologies_present']['SSO']:
                                if gto['ontologies_present']['SSO'][term] != 'Unknown':
                                    term = gto['ontologies_present']['SSO'][term]

                        # convert terms to rxns
                        if term in translations[ontology_type]:
                            rxn = translations[ontology_type][term]
                        else:
                            if oe in summary["orphan_terms"]:
                                summary["orphan_terms"][oe].append(term)
                                summary["orphan_terms"][oe] = list(set(summary["orphan_terms"][oe]))
                            else:
                                summary["orphan_terms"][oe] = [term]

                        # terms
                        if term in summary["genes"][gene_id]['terms']:
                            summary["genes"][gene_id]['terms'][term].append(oe)
                        else:
                            summary["genes"][gene_id]['terms'][term] = [oe]

                        if term in summary['terms']:
                            summary['terms'][term].append(oe)
                            summary['terms'][term] = list(set(summary['terms'][term]))
                        else:
                            summary['terms'][term] = [oe]

                        # rxns
                        if rxn != "none":
                            if rxn in summary["genes"][gene_id]['rxns']:
                                summary["genes"][gene_id]['rxns'][rxn].append(oe)
                            else:
                                summary["genes"][gene_id]['rxns'][rxn] = [oe]

                            if rxn in summary['rxns']:
                                summary['rxns'][rxn].append(oe)
                                summary['rxns'][rxn] = list(set(summary['rxns'][rxn]))
                            else:
                                summary['rxns'][rxn] = [oe]

    return summary


gto = json.loads(open("PDIF272563.7.json", "r").read())
translations = json.loads(open("translations.json", "r").read())

summary = summarize(gto, translations)

with open('summary.json', 'w') as outfile:
    json.dump(summary, outfile, indent=2)


# format

html_summary_report = {}

for oe in summary['ontology_events']:
    html_summary_report[oe] = {"gene": [], "term": [], "rxn": []}

for gene in summary["genes"]:
    for term in summary["genes"][gene]['terms']:
        for oe in summary["genes"][gene]['terms'][term]:
            html_summary_report[oe]['gene'].append(gene)
            html_summary_report[oe]['term'].append(term)

            html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))
            html_summary_report[oe]['term'] = list(set(html_summary_report[oe]['term']))

    for rxn in summary["genes"][gene]['rxns']:
        for oe in summary["genes"][gene]['rxns'][rxn]:
            html_summary_report[oe]['rxn'].append(rxn)
            html_summary_report[oe]['gene'].append(gene)

            html_summary_report[oe]['rxn'] = list(set(html_summary_report[oe]['rxn']))
            html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))

with open('html_summary.json', 'w') as outfile2:
    json.dump(html_summary_report, outfile2, indent=2)
