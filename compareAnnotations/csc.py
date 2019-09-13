import json
import pandas as pd

from bokeh.plotting import figure, output_file
from bokeh.io import show
from bokeh.palettes import viridis
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.transform import factor_cmap


def read_summary(summary_file):
    return json.loads(open(summary_file, "r").read())


def get_descriptions(summary):
    d = {}
    for o in summary["ontology_events"]:
        d[o] = summary["ontology_events"][o].get(
            'description', summary["ontology_events"][o]['method']) + '_' + str(o)
    return(d)


def longest_set(s, w):
    s = s.copy()
    for event in s:
        for winner in w:
            s[event] = s[event] - w[winner]

    # https://stackoverflow.com/a/21839239
    max_key, max_value = max(s.items(), key=lambda x: len(x[1]))
    return(max_key)


def plot_csc2(summary, summary_type="rxn"):
    descriptions = get_descriptions(summary)
    events = sorted(summary['ontology_events'].keys())
    rxns = summary[summary_type]

    # convert to sets
    rxns_in_events = dict((int(el), set()) for el in events)
    for rxn in rxns:
        for event in rxns[rxn]:
            rxns_in_events[event].add(rxn)

    winning_sets = {}
    winning_order = []
    baseline = 0
    df = pd.DataFrame(columns=["E", "C", "T", "L", "R"])
    # E=event, C=comparison, T=total, L=left, R=right

    for _ in range(len(rxns_in_events)):

        current_right = baseline
        current_left = baseline

        # get current winner
        longest_set_key = longest_set(rxns_in_events, winning_sets)

        # compare current winner to all past winners
        current = rxns_in_events[longest_set_key]
        for past_winner in winning_order:
            overlap = len(winning_sets[past_winner] & current)
            current_left -= overlap
            row = [descriptions[str(longest_set_key)],  # E
                   descriptions[str(past_winner)],  # C
                   overlap,  # T
                   current_left,  # L
                   current_left + overlap]  # R
            df.loc[len(df)] = row
            current = current - winning_sets[past_winner]

        # process current winner
        row = [descriptions[str(longest_set_key)],  # E
               descriptions[str(longest_set_key)],  # C
               len(current),  # T
               current_right,  # L
               current_right + len(current)]  # R

        df.loc[len(df)] = row  # add to df
        baseline += len(current)

        # move current winner to past winners
        winning_sets[longest_set_key] = rxns_in_events[longest_set_key]
        winning_order.append(longest_set_key)
        rxns_in_events[longest_set_key] = set()

    source = ColumnDataSource(df)

    type1_colormap = factor_cmap('E', palette=viridis(
        len(df.E.unique())), factors=df.E.unique())
    type2_colormap = factor_cmap('C', palette=viridis(
        len(df.C.unique())), factors=df.C.unique())

    p = figure(y_range=df.E.unique().tolist()[::-1],  # .tolist()[::-1] reverses the list.
               plot_height=300,
               plot_width=1000,
               title="Annotation events ranked by \'" + str(summary_type) + "\' contribution",
               tools="wheel_zoom,box_zoom,reset,save")

    p.hbar(y='E',
           height=0.9,
           left='L',
           right='R',
           source=source,
           fill_color=type2_colormap,
           line_color="black")

    p.add_tools(HoverTool(tooltips=[("Total", "@T"), ("Comparison", "@C")]))
    p.title.text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"
    return p


# summary = read_summary("PT19DW.5.json")
summary = read_summary("PT19DW.7.json")

csc_rxns2 = plot_csc2(summary, "rxns")
#csc_rxns = plot_csc(summary, "rxn")
#
output_file("csc_rxns.html", title="CSC")
show(csc_rxns2)
