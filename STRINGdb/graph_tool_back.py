from typing import Union,Literal,Optional,Iterable
from pathlib import Path
import gzip
import csv
import pandas as pd
import matplotlib.colors as mpc
import matplotlib.cm as cm
import graph_tool.all as gt
import igraph as ig
from . import __version__, prep
from .prep import Identifier, DbFile




class Graph:
    def __init__(self, idents: Identifier):
        self.species = idents.species
        self.ids = idents.ids

    def subnetwork(self):
        graphfile = f'{self.species}.{DbFile().links}.txt.gz'
        if not Path.joinpath(prep.cache_path,graphfile).is_file():
            DbFile().download(self.species,'links')
        with gzip.open(Path.joinpath(prep.cache_path,graphfile), 'rt') as f:
            graph = gt.load_graph_from_csv(f, directed=True,
                    skip_first=True,csv_options={"delimiter": " "})
        subg = gt.GraphView(graph, vfilt=self.ids)

