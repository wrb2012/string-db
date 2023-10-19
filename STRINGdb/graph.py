from typing import List,Dict,Union,Literal,Optional,Iterable
from pathlib import Path
import gzip
import pandas as pd
import matplotlib.colors as mpc
import matplotlib.cm as cm
import igraph as ig
from . import prep
from .prep import Identifier, DbFile




class Graph:
    def __init__(self, idents: Identifier):
        self.species = idents.species
        self.ids = idents.ids

    def subnetwork(self, score_source: Literal ='combined_score') -> pd.DataFrame:
        '''return: (string_ids, interactions)'''
        graphfile = f'{self.species}.{DbFile().links}.txt.gz'
        if not Path.joinpath(prep.cache_path,graphfile).is_file():
            DbFile().download(self.species,'links')
        ref = pd.read_csv(Path.joinpath(prep.cache_path,graphfile), sep=' ')
        refgraph = ig.Graph.TupleList(ref[['protein1','protein2', score_source]].itertuples(index=False),
                directed=True, edge_attrs=score_source)
        subgraph: ig.Graph = refgraph.induced_subgraph(self.ids)
        interactions = subgraph.to_dict_list(use_vids=False)[1]
        return pd.DataFrame(interactions)
    
    def plot(self, type = 'stress'):
        pass

    def neighbors(self):
        pass