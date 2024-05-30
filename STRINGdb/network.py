from io import StringIO
from typing import Dict, Iterable, Literal, Optional, Union

import matplotlib.cm as cm
import matplotlib.colors as mpc
import numpy as np
import pandas as pd

from . import __version__, prep
from .prep import Identifier


class Image:
    """not recommend plotting over 150 genes/proteins
    use .params() first"""

    def __init__(self, idents: Identifier):
        self.ids = idents.ids
        self.sig = idents.sig
        self.data = {
            "species": idents.species,
            "identifiers": idents(),
            "caller_identity": __package__,
        }

    def params(self,
        thres: int = 400,
        network_type: str = "functional",
        edge_style: Literal["confidence", "evidence"] = "evidence",
        hide_disconnect=0,
        flat_node=0,
        center_label=0,
        label_size=12,
    ):
        self.data["required_score"] = thres
        if edge_style == "confidence":
            self.data["network_flavor"] = edge_style
        self.data["flat_node_design"] = flat_node
        if label_size != 12:
            self.data["custom_label_font_size"] = label_size

    def add_holo(self, logFoldChange: str = "logFC") -> pd.DataFrame:
        anno = pd.DataFrame(self.sig)
        pos = anno[anno[logFoldChange] > 0].reset_index()
        norm = np.exp(pos[logFoldChange]) / np.exp(pos[logFoldChange].max())
        color_pos = cm.ScalarMappable(cmap=cm.Reds).to_rgba(norm)
        color_pos_hex = pd.Series([mpc.to_hex(c) for c in color_pos], name='color')
        pos_df = pos.join(color_pos_hex)  #if no .reset_index(), join will fail
        neg = anno[anno[logFoldChange] <= 0].reset_index()
        norm = np.exp(neg[logFoldChange]) / np.exp(neg[logFoldChange].max())
        color_neg = cm.ScalarMappable(cmap=cm.Blues).to_rgba(norm)
        color_neg_hex = pd.Series([mpc.to_hex(c) for c in color_neg], name='color')
        neg_df = neg.join(color_neg_hex)
        anno_colors = pd.concat([pos_df, neg_df])

        self.data["colors"] = anno_colors["color"].str.cat(sep="\r")
        res = prep.client.post("/cgi/webservices/post_payload.pl", data=self.data)

        self.data["internal_payload_id"] = res.text
        return anno_colors

    def plot(self,
        img: Literal["svg", "png"] = "svg", *, save: str = ""
    ) -> Union[bytes, str, None]:
        """params:
        img : svg or png

        return:

        """
        _format = "highres_image" if img == "png" else "svg"
        res = prep.client.post("/api/" + _format + "/network", data=self.data)
        image = res.text if img == "svg" else res.content
        if save:
            if img == "svg":
                with open(save, "w") as f:
                    f.write(image)
            elif img == "png":
                with open(save, "wb") as f:
                    f.write(image)
        else:
            return image


class Enrichment:
    def __init__(self, idents: Identifier):
        self.idents = idents
        self.sig = idents.sig
        self.data = {
            "species": idents.species,
            "identifiers": idents(),
            "caller_identity": __package__ + "/" + __version__,
        }

    def params(
        self,
        thres: int = 400,
        background: Iterable = [],
    ):
        self.data["required_score"] = thres
        self.data["background"] = "\r".join(background)

    def interaction(self) -> pd.DataFrame:
        res = prep.client.post("/api/tsv/network", data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep="\t")

    def all_partner(self, limit=None) -> pd.DataFrame:
        if limit:
            self.data["limit"] = limit
        res = prep.client.post("/api/tsv/interaction_partners", data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep="\t")

    def similarity(self) -> pd.DataFrame:
        res = prep.client.post("/api/tsv/homology", data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep="\t")

    def ortholog(self, other_species: Iterable) -> pd.DataFrame:
        if not other_species:
            other_species = [self.idents.species]
        data = {"species_b": "\r".join(other_species)}.update(self.data)
        res = prep.client.post("/api/tsv/homology_best", data=data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep="\t")

    def functional(self,
        category: Literal["All, Process, Component, Function, Keyword, KEGG, RCTM, Pfam, SMART, InterPro"] = "All",
    ) -> pd.DataFrame:
        res = prep.client.post("/api/tsv/enrichment", data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep="\t")

    def functional_annotation(self, allow_pubmed: bool = False) -> pd.DataFrame:
        if allow_pubmed:
            self.data["allow_pubmed"] = "1"
        res = prep.client.post("/api/tsv/functional_annotation", data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep="\t")

    def ppi(self, thres: Optional[int] = None) -> Dict:
        if thres:
            self.data["required_score"] = thres
        # if 'background_string_identifiers' not in self.data:
        res = prep.client.post("/api/json/ppi_enrichment", data=self.data)
        return res.json()[0]
