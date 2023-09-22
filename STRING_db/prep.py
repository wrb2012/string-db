
from enum import Enum
from pathlib import Path
from typing import Dict,List,Set,Union
import httpx
import pandas as pd
from .version import api_domain, STRING_VER, STATIC_FILES

client = httpx.Client(base_url=api_domain+'/api',
        timeout=20,
        transport=httpx.HTTPTransport(retries=2)
)


class Identifier:

    def __init__(self, taxon: int, idents: Union[List,Set] = {}):
        self.species = taxon
        self.idents = idents
    
    def __call__(self) -> str:
        if len(self.idents) > 1:
            return '\r'.join(self.idents)  # for POST method
        else:
            return self.idents[0]

class Table(str, Enum):
    '''database resources'''
    info = 'protein.info.v'+STRING_VER
    alias = 'protein.aliases.v'+STRING_VER
    links = 'protein.links.detail.v'+STRING_VER
    links_full = 'protein.links.full.v'+STRING_VER
    physical = 'protein.physical.links.detail.v'+STRING_VER
    physical_full = 'protein.physical.links.detail.v'+STRING_VER

    sequence = 'protein.sequence.v'+STRING_VER
    homology = 'protein.homology.v'+STRING_VER
    enrich = 'protein.enrichment.v'+STRING_VER
    cluster = 'clusters.info.v'+STRING_VER
    cluster_tree = 'clusters.tree.v'+STRING_VER

default_cache = Path.joinpath(Path.home(),'.cache', __package__)

def download_table(taxon: int, table: str, save_dir: Path = default_cache):
    '''download large database table'''
    Path.mkdir(Path(save_dir), exist_ok=True)
    
    url_lastest = f'{STATIC_FILES}{Table[table].value}/{taxon}.{Table[table].value}.txt.gz'
    with open(f'{save_dir}/{taxon}.{Table[table].value}.txt.gz', 'wb') as f:
        with httpx.stream("GET", url_lastest) as res:
            for chunk in res.iter_raw():
                #file = gzip.decompress(res.content)
                f.write(chunk)


class IDConvert:
    
    def __init__(self, taxon: int = 9606):
        self.species = taxon

    def map_id(self, idents: str) -> List[Dict]:
        '''returned fields:
        queryIndex	position of the protein in your input (starting from position 0)
        stringId	STRING identifier
        ncbiTaxonId	NCBI taxon identifier
        taxonName	species name
        preferredName	common protein name
        annotation	protein annotation
        '''
        res = client.post('/json/get_string_ids',
                data={'identifiers': idents, 'species': self.species})
        return res.json()  ## usually you just need first ids
    
    def map_id_local(source: Union[Set,List], ref_dir = default_cache):
        protein_aliases: pd.DataFrame = pd.read_csv(ref_dir, sep='\t')
        ## warning : may loss genes/proteins
        _to_stringid = protein_aliases.drop_duplicates(subset=['alias','#string_protein_id'])
        string_ids = _to_stringid[_to_stringid['alias'].isin(source)]\
                .groupby(['alias'])['#string_protein_id'].first()
        print(f'Conversion rate is: {len(string_ids) / len(source)} !')
        return string_ids.values
