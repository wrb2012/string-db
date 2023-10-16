from typing import Union,Set,List,Literal
import numpy as np
import pandas as pd
import matplotlib.colors as mpc
import matplotlib.cm as cm
from . import __version__
from .prep import client, Identifier

class Image:
    '''not recommend plotting over 150 genes/proteins
    use .params() first'''
    def __init__(self,idents: Identifier):
        self.idents = idents
        self.sig = idents.sig
        self.data = {"species": idents.species,
                     'caller_identity': __package__+'/'+__version__}

    def params(self, de_color: dict = None, thres = 400, network_type: str = 'functional',
               edge_style = "confidence", hide_disconnect = 0, flat_node = 0,
               center_lable = 0, lable_size = 12):
        self.data['identifiers'] = self.idents()
        self.data['required_score'] = thres
        #self.data['add_color_nodes'] = de_color
        
    def visual(self, logFoldChange: str='logFC') -> pd.DataFrame:
        
        anno = pd.DataFrame(self.sig[logFoldChange])
        anno = anno.reset_index(drop=True)   ### for slice input
        norm = mpc.Normalize(vmin=anno.min(), vmax=anno.max())
        color_array = cm.ScalarMappable(cmap=cm.bwr, norm=norm).to_rgba(anno)
        color_hex = pd.Series([mpc.to_hex(c, keep_alpha=True).upper()
                for c in color_array], name='color')
        anno_colors = anno.join(color_hex)

        #df2 = pd.concat([pos_df, neg_df])
        #df2.reindex(anno.index)
        
        data = {'identifiers': '\r'.join(self.idents.ids)}
        data['color'] = anno_colors['color'].str.cat(sep='\r')
        ## 这里服务器返回所需时间较长
        res = client.post('/cgi/webservices/post_payload.pl',
            data=data)

        self.data['internal_payload_id'] = res.text
        return anno_colors

    def network(self, img: Literal['svg','png'] = 'svg', *, save: Union[bool,str]=True
    ) -> Union[bytes,str]:
        '''params:
        img : svg or png
        
        return:
        
        '''
        _format = 'highres_image' if img=='png' else 'svg'
        res = client.post('/api/'+_format+'/network', data=self.data)
        
        if img == 'svg':
            return res.text
        elif img == 'png':
            return res.content

