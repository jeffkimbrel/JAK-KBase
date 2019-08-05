import json

from bokeh.plotting import figure, output_file
from bokeh.io import show
from bokeh.palettes import inferno
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.transform import factor_cmap
from bokeh.models import HoverTool
# from bokeh.io import export_svgs


def read_summary(summary_file):
    return json.loads(open(summary_file, "r").read())


def get_descriptions(summary):
    d = {}
    for o in summary["ontology_events"]:
        print(o)
        d[o] = summary["ontology_events"][o].get(
            'description', summary["ontology_events"][o]['method']) + '_' + str(o)
    return(d)


def plot_totals(summary):
    descriptions = get_descriptions(summary)
    totals = {}
    for event in summary['ontology_events'].keys():
        totals[str(event)] = {'genes': [],
                              'rxns': [],
                              'terms': []}

    # genes
    for gene in summary['genes']:
        for term in summary['genes'][gene]['terms']:
            for event in summary['genes'][gene]['terms'][term]:
                totals[str(event)]['genes'].append(gene)

    # terms
    for term in summary['terms']:
        for event in summary['terms'][term]:
            totals[str(event)]['terms'].append(term)

    # rxns
    for rxn in summary['rxns']:
        for event in summary['rxns'][rxn]:
            totals[str(event)]['rxns'].append(rxn)

    # sums
    events = []
    types = ['genes', 'terms', 'rxns']

    gene_counts = []
    rxn_counts = []
    term_counts = []

    for event in totals:
        events.append(descriptions[event])
        gene_counts.append(len(set(totals[event]['genes'])))
        rxn_counts.append(len(set(totals[event]['rxns'])))
        term_counts.append(len(set(totals[event]['terms'])))

    data = {'events': events,
            'genes': gene_counts,
            'terms': term_counts,
            'rxns': rxn_counts
            }

    x = [(event, type) for event in events for type in types]

    counts = sum(zip(data['genes'], data['terms'], data['rxns']), ())
    source = ColumnDataSource(data=dict(x=x, counts=counts))

    p = figure(y_range=FactorRange(*x),
               plot_height=400,
               plot_width=1000,
               title="Unique Counts per Annotation Event",
               tools="wheel_zoom,box_zoom,reset,save")

    p.hbar(y='x',
           right='counts',
           height=0.9,
           source=source,
           line_color="black",
           fill_color=factor_cmap('x',
                                  palette=inferno(len(types)),
                                  factors=types,
                                  start=1,
                                  end=2))

    p.x_range.start = 0
    p.y_range.range_padding = 0.1
    p.yaxis.major_label_orientation = "horizontal"
    p.yaxis.subgroup_label_orientation = "horizontal"
    p.yaxis.group_label_orientation = "horizontal"
    p.ygrid.grid_line_color = None
    p.title.text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"
    p.yaxis.group_text_font_size = "12pt"
    p.add_tools(HoverTool(tooltips=[("Type", "@x"), ("Count", "@counts")]))

    return(p)


#summary = read_summary("PT19DW.5.json")
summary = read_summary("summary.json")

output_file("totals.html", title="Totals")
totals = plot_totals(summary)

show(totals)
