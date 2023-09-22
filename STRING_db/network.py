from typing import Union,Set,List
import pandas as pd
from . import __version__
from .prep import client, Identifier

class Image:
    
    def __init__(self,idents: Identifier):
        self.idents = idents
        self.data = {"species": idents.species}

    def params(self, de_color: dict = None, thres = 400, network_type: str = 'functional'):
        self.data['caller_identity'] = __package__+'/'+__version__
        self.data['identifiers'] = self.idents()
        self.data['required_score'] = thres
        self.data['add_color_nodes'] = de_color
        
    def visual():
        pass
         
    def network(self, img: str = 'svg') -> Union[bytes,str]:
        '''params:
        img : svg or png
        
        return:
        
        '''
        _format = 'highres_image' if img=='png' else 'svg'
        res = client.post(_format+'/network', data=self.data)
        
        if img == 'svg':
            return res.text
        elif img == 'png':
            return res.content

