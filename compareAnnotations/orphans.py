import json

from bokeh.plotting import figure, output_file
from bokeh.io import show
from bokeh.palettes import viridis
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
# from bokeh.io import export_svgs


def read_summary(summary_file):
    return json.loads(open(summary_file, "r").read())


def get_descriptions(summary):
    d = {}
    for o in summary["ontology_events"]:
        d[o] = summary["ontology_events"][o].get(
            'description', summary["ontology_events"][o]['method']) + '_' + str(o)
    return(d)


def plot_orphans(summary):
    descriptions = get_descriptions(summary)
    events = sorted(summary['ontology_events'].keys())
    orphan_terms = summary['orphan_terms']
    orphans = []

    for event in events:
        if event in orphan_terms:
            orphans.append(len(set(orphan_terms[event])))
        else:
            orphans.append(0)

    events_mapped = []
    for event in events:
        events_mapped.append(descriptions[event])

    source = ColumnDataSource(data=dict(events=events_mapped,
                                        orphans=orphans,
                                        color=viridis(len(events))))

    p = figure(y_range=events_mapped,
               plot_height=300,
               plot_width=1000,
               title="Terms without a modelSEED Reaction",
               tools="wheel_zoom,box_zoom,reset,save")

    p.hbar(y='events',
           height=0.7,
           left=0,
           right='orphans',
           color='color',
           source=source,
           line_color="black")

    p.add_tools(HoverTool(tooltips=[("Orphan Terms", "@orphans")]))
    p.title.text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"
    return(p)


#summary = read_summary("PT19DW.5.json")
summary = read_summary("summary.json")

orphans = plot_orphans(summary)

output_file("orphans.html", title="Orphans")
show(orphans)
