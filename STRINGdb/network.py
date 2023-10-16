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
        
    def visual(self, logFoldChange: str='logFC') -> bytes:
        
        anno = pd.DataFrame(self.sig)
        pval05_pos = anno[anno[logFoldChange] > 0].reset_index(drop=True)
        norm = np.exp(pval05_pos[logFoldChange]) / np.exp(pval05_pos[logFoldChange].max())
        color_pos = cm.ScalarMappable(cmap=cm.Reds).to_rgba(norm)
        color_pos_hex = pd.Series([mpc.to_hex(c, keep_alpha=True) for c in color_pos],
                name='color')
        pos_df = pval05_pos.join(color_pos_hex)  #if no .reset_index(), join will fail

        pval05_neg = anno[anno[logFoldChange] <= 0].reset_index(drop=True)
        norm = np.exp(pval05_neg[logFoldChange].abs()) / np.exp(-pval05_neg[logFoldChange].min())
        color_neg = cm.ScalarMappable(cmap=cm.Greens).to_rgba(norm)
        color_neg_hex = pd.Series([mpc.to_hex(c, keep_alpha=True) for c in color_neg],
                name='color')
        neg_df = pval05_neg.join(color_neg_hex)

        color_df = pd.concat([pos_df, neg_df])
        #color_df.reindex(anno.index)
        
        data = {'identifiers': '\r'.join(self.idents.ids)}
        data['color'] = color_df['color'].str.cat(sep='\r')
        ## 这里服务器返回所需时间较长
        res = client.post('/cgi/webservices/post_payload.pl',
            data=data)

        self.data['internal_payload_id'] = res.text
        return color_df

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

