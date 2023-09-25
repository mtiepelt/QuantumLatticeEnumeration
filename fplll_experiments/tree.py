import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib.pyplot as plt
import os


def extract_graph(node):
    graph = []
    node_lbl = node['n'] + "|" + str(node['v'])
    for child in node['c']:
        child_lbl = child['n'] + "|" + str(child['v'])
        graph.append((node_lbl, child_lbl))
        graph += extract_graph(child)
    return graph


def draw_graph(tree, title=None, filename="tree.png", node_text_size=12, text_font='sans-serif'):

    # extract graph from fplll tree
    graph = extract_graph(tree)

    # create networkx graph
    G=nx.DiGraph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    if isinstance(title, str):
        plt.title(title)

    # graph_pos = graphviz_layout(G, prog='dot', args='-Gnodesep=3')
    # nx.draw(G, graph_pos, with_labels=False, arrows=True)

    # nx.draw_networkx_labels(G, graph_pos, font_size=node_text_size,
    #                         font_family=text_font)

    # # save plot
    # plt.savefig(filename)

    # # write dot file to use with graphviz
    # # run "dot -Tpng test.dot >test.png"
    write_dot(G, f'{filename}.dot')

    os.system(f"dot -Tpng {filename}.dot >{filename}")
    os.system(f"rm -rf {filename}.dot")


def subtree_to_vals(node, lbl="0"):
    vals = [(lbl, node['v'])]
    cidx = 0
    for child in node['c']:
        vals += subtree_to_vals(child, lbl=f"{lbl}-{cidx}")
        cidx += 1
    return vals


def num_nodes(node):
    return len(subtree_to_vals(node))


def extract_distribution(vals):
    distr = {}
    for key, val in vals:
        lvl = key.count('-')
        if lvl not in distr:
            distr[lvl] = []
        distr[lvl].append(val)
    return distr

def subtree_to_nodes_per_level(node):
    vals = subtree_to_vals(node, "0")
    levels = {}
    for lbl, val in vals:
        lvl = lbl.count('-')
        if lvl not in levels:
            levels[lvl] = 0
        levels[lvl] += 1
    return levels

def main():
    import argparse
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument('tree', type=str, help='tree json file')
    args = parser.parse_args()

    with open(args.tree) as f:
        bkz_s = "".join(f.readlines())
    bkz = json.loads(bkz_s)
    # d = extract_distribution(subtree_to_vals(tree))
    # print(d)

    for tour in range(len(bkz)):
        for index in bkz[tour]:
            print((tour, index))
            tree = bkz[tour][index]
            draw_graph(tree, title=f"tour {tour}, index {index}", filename=f"{tour}-{index}.png")


if __name__ == '__main__':
    main()



