from typing import List,Dict,Union,Literal,Optional,Iterable
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import igraph as ig
from . import prep
from .prep import Identifier, DbFile



class Graph:
    def __init__(self, idents: Identifier):
        self.species = idents.species
        self.ids = idents.ids

    def subnetwork(self, score_source: Literal ='combined_score', /, thres = None,
        extra_label: Optional[pd.Series]=[]) -> pd.DataFrame:
        '''return: interactions(subset network dataframe)'''
        graphfile = f'{self.species}.{DbFile().links}.txt.gz'
        if not Path.joinpath(prep.cache_path,graphfile).is_file():
            DbFile().download(self.species,'links')
        ref = pd.read_csv(Path.joinpath(prep.cache_path,graphfile), sep=' ')
        refgraph = ig.Graph.TupleList(ref[['protein1','protein2', score_source]].itertuples(index=False),
                directed=True, edge_attrs='score')
        self.subgraph: ig.Graph = refgraph.induced_subgraph(self.ids)
        if thres:
            edges = self.subgraph.es(score_gt=thres)
            self.subgraph = self.subgraph.subgraph_edges(edges)
        if len(extra_label) > 0:
            self.subgraph.vs['label'] = [extra_label[extra_label == s].index[0]
                    for s in self.subgraph.vs['name']]
        interactions = self.subgraph.to_dict_list(use_vids=False)[1]
        return pd.DataFrame(interactions)
    
    def plot(self, layout = 'kk', *, label: bool=True) -> plt.Figure:
        '''ref: <https://python.igraph.org/en/stable/tutorial.html#layouts-and-plotting>
        '''
        components = self.subgraph.connected_components(mode='weak')
        fig, ax = plt.subplots()
        ig.plot(components, ax, (0,0,1200,1000),
            layout=self.subgraph.layout(layout),
            palette=ig.RainbowPalette(),
            vertex_size=[0.1+ s/50 for s in self.subgraph.degree()],
            vertex_color=list(map(int, ig.rescale(components.membership, (0, 200), clamp=True))),
            vertex_label=self.subgraph.vs['label'],
            edge_width=[0.1+ s/1000 for s in self.subgraph.es['score']],
        )
        return fig.figure
    def neighbors(self):
        pass