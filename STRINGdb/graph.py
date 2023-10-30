from typing import List,Dict,Union,Literal,Optional,Iterable
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import igraph as ig
from . import prep
from .prep import Identifier, DbFile



class Graph:
    def __init__(self, idents: Identifier, score_source: Literal ='combined_score'):
        self.species = idents.species
        self.ids = idents.ids
        graphfile = DbFile().download(idents.species,'links')
        ref = pd.read_csv(graphfile, sep=' ')
        self.refgraph = ig.Graph.TupleList(ref[['protein1','protein2', score_source]].itertuples(index=False),
                directed=True, edge_attrs='score')
    
    def subnetwork(self, /, thres = None, label: Optional[pd.Series]=[]) -> pd.DataFrame:
        '''return: interactions(subset network dataframe)'''
        
        self.subgraph: ig.Graph = self.refgraph.induced_subgraph(self.ids)
        if thres:
            edges = self.subgraph.es(score_gt=thres)
            self.subgraph = self.subgraph.subgraph_edges(edges)
        if len(label) > 0:
            self.subgraph.vs['label'] = [label[label == s].index[0]
                    for s in self.subgraph.vs['name']]
        interactions = self.subgraph.to_dict_list(use_vids=False)[1]
        self.interactions = pd.DataFrame(interactions)
        return self.interactions
    
    def plot(self, graph='all', /, layout = 'auto', *, vertex_label: bool=True) -> plt.Figure:
        '''ref: <https://python.igraph.org/en/stable/tutorial.html#layouts-and-plotting>
        '''
        deduped = self.subgraph.copy() if graph == 'all' else self.nonHomologous.copy()
        deduped.to_undirected('collapse', combine_edges='first')
        components = deduped.connected_components(mode='weak')
        fig, ax = plt.subplots(dpi=300) #figsize=(10, 10)
        ig.plot(components, ax,
            layout=self.subgraph.layout(layout),
            palette=ig.RainbowPalette(),
            vertex_size=[4+2*s for s in self.subgraph.degree()],
            vertex_color=list(map(int, ig.rescale(components.membership, (0, 200), clamp=True))),
            vertex_label=self.subgraph.vs['label'],
            vertex_label_dist=0.5,
            vertex_label_size=5,
            edge_width=[0.125+ w/1000 for w in self.subgraph.es['score']],
        )
        ax.set_title('Clutered Network')
        return fig

    def get_neighbors(self) -> List[str]:
        neighbor_vertices = []
        for id in self.ids:
            neighbor_vertices += self.refgraph.neighbors(id)
        neighbor_ids = self.refgraph.vs(neighbor_vertices)
        return list(set(neighbor_ids))

    def clustreing(self, algorithm=None):
        if not hasattr(self, 'subgraph'):
            self.subnetwork()
        if algorithm=="fastgreedy":
            communities = self.subgraph.community_fastgreedy()
        elif algorithm=="fastgreedy":
            communities = self.subgraph.community_walktrap()
        elif algorithm=="leiden":
            communities = self.subgraph.community_leiden()
        else:
            communities = self.subgraph.community_edge_betweenness()
        communities = communities.as_clustering()
        return communities

    def homology_ref(self, min_bitscore: int=60) -> ig.Graph:
        homology_file = DbFile().download(self.species,'homology')
        homology = pd.read_csv(Path.joinpath(prep.cache_path,homology_file), sep='\t')
        if min_bitscore:
            homology = homology[homology['bitscore']>min_bitscore]
        homology_graph = ig.Graph.TupleList(homology.iloc[:, [0,1,3]].itertuples(index=False),
                directed=True, edge_attrs='bitscore')
        homoref = homology_graph.induced_subgraph(self.ids)
        return homoref
    
    def homologous_interactions(self, interactions_graph: ig.Graph=None, thres: int=60):
        if not interactions_graph:
            interactions_graph = self.subgraph
        self.homologous = interactions_graph.intersection(self.homology_ref())
        homologous_interactions = self.homologous.to_dict_list(use_vids=False)[1]
        return pd.DataFrame(homologous_interactions)

    def remove_homologous_interactions(self, /,
        thres: int= 60):
        hg = self.homologous_interactions(thres=thres) 
        self.nonHomologous = self.subgraph.difference(self.homology_ref()) #bug when using difference
        print(f'{len(hg)} interaction(s) removed')
        tmp = pd.concat([self.interactions, hg[['source','target', 'score']]])
        return tmp.drop_duplicates(keep=False)

