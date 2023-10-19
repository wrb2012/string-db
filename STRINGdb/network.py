from typing import Union,Literal,Optional,Iterable
from io import StringIO
import pandas as pd
import matplotlib.colors as mpc
import matplotlib.cm as cm
from . import __version__, prep
from .prep import Identifier


class Image:
    '''not recommend plotting over 150 genes/proteins
    use .params() first'''
    def __init__(self,idents: Identifier):
        self.idents = idents
        self.sig = idents.sig
        self.data = {"species": idents.species,
                     'caller_identity': __package__+'/'+__version__}

    def params(self, thres = 400, network_type: str = 'functional',
               edge_style: Literal['confidence','evidence']='evidence', 
               hide_disconnect = 0, flat_node = 0,
               center_label = 0, label_size = 12):
        self.data['identifiers'] = self.idents()
        self.data['required_score'] = thres
        if edge_style == 'confidence':
            self.data['network_flavor'] = edge_style
        self.data['flat_node_design'] = flat_node
        if label_size != 12:
            self.data['custom_label_font_size'] = label_size
    def add_holo(self, logFoldChange: str='logFC') -> pd.DataFrame:
        
        anno = pd.DataFrame(self.sig[logFoldChange])
        anno = anno.reset_index(drop=True)   ### for slice input
        norm = mpc.Normalize(vmin=anno.min(), vmax=anno.max())
        color_array = cm.ScalarMappable(cmap=cm.bwr, norm=norm).to_rgba(anno)
        color_hex = pd.Series([mpc.to_hex(c, keep_alpha=True).upper()
                for c in color_array], name='color')
        anno_colors = anno.join(color_hex)

        #df2 = pd.concat([pos_df, neg_df])
        #df2.reindex(anno.index)
        
        data = {'identifiers': ' '.join(self.idents.ids)}
        data['colors'] = anno_colors['color'].str.cat(sep=' ')
        res = prep.client.post('/cgi/webservices/post_payload.pl',
            data=data)

        self.data['internal_payload_id'] = res.text
        return anno_colors

    def plot(self, img: Literal['svg','png']='svg', *, save: Union[bool,str]=True
    ) -> Union[bytes,str]:
        '''params:
        img : svg or png
        
        return:
        
        '''
        _format = 'highres_image' if img=='png' else 'svg'
        res = prep.client.post('/api/'+_format+'/network', data=self.data)
        
        if img == 'svg':
            return res.text
        elif img == 'png':
            return res.content



class Enrichment:
    def __init__(self,idents: Identifier):
        self.idents = idents
        self.sig = idents.sig
        self.data = {"species": idents.species,
                     'identifiers': self.idents(),
                     'caller_identity': __package__+'/'+__version__}

    def params(self, thres: int= 400, background: Iterable=[],):
        self.data['required_score'] = thres
        self.data['background'] = '\r'.join(background)
    
    def functional(self, category: Literal['All, Process, Component, Function, Keyword, KEGG, RCTM, Pfam, SMART, InterPro']='All'):
        res = prep.client.post('/api/tsv/enrichment', data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep='\t')
    
    def functional_annotation(self, allow_pubmed:bool=False):
        if allow_pubmed:
            self.data['allow_pubmed'] = '1'
        res = prep.client.post('/api/tsv/functional_annotation', data=self.data)
        tsv = StringIO(res.text)
        return pd.read_csv(tsv, sep='\t')
    
    def ppi(self , thres: Optional[int]=None):
        if thres:
            self.data['required_score'] = thres
        #if 'background_string_identifiers' not in self.data:
        res = prep.client.post('/api/json/ppi_enrichment', data=self.data)
        return res.json()
