
from pathlib import Path
from io import StringIO
import sys
from typing import Dict,List,Set,Union,Iterable,Optional
import httpx
import pandas as pd
from . import __version__
from .version import choose_version, STATIC_ASSET

default_cache = Path.joinpath(Path.home(),'.cache', __package__)

def init(version: float=12.0, *, ref_dir=default_cache):
    domain = choose_version(version)
    global client
    client = httpx.Client(base_url=domain,
            headers={'User-Agent': __package__+'/'+__version__},
            timeout=20,
            transport=httpx.HTTPTransport(retries=3)
    )
    global STRING_VER
    STRING_VER = version
    global cache_path
    cache_path = Path(ref_dir)
    

class Identifier:

    def __init__(self, taxon: int, idents: Iterable = [], *,
            sig = None):
        self.species = taxon
        self.ids = idents
        self.sig = sig

    def __call__(self) -> str:
        if len(self.ids) > 1:
            return '\r'.join(self.ids)  # for POST method
        else:
            return self.ids[0]


class DbFile:
    '''database resources'''
    def __init__(self, string_ver: Optional[str]=None):
        string_ver = str(string_ver) if string_ver else str(STRING_VER)
        self.info = 'protein.info.v'+string_ver
        self.alias = 'protein.aliases.v'+string_ver
        self.links = 'protein.links.detailed.v'+string_ver
        self.links_full = 'protein.links.full.v'+string_ver
        self.physical = 'protein.physical.links.detailed.v'+string_ver
        self.physical_full = 'protein.physical.links.full.v'+string_ver

        self.sequence = 'protein.sequence.v'+string_ver
        self.homology = 'protein.homology.v'+string_ver
        self.enrich = 'protein.enrichment.v'+string_ver
        self.cluster = 'clusters.info.v'+string_ver
        self.cluster_tree = 'clusters.tree.v'+string_ver

    def download(self, taxon: int, table: str, *, save_dir: Path = default_cache) -> Path:
        '''download large database table'''
        if cache_path:
            save_dir = cache_path
        Path.mkdir(Path(save_dir), exist_ok=True)
        file = Path.joinpath(save_dir, f'{taxon}.{getattr(self,table)}.txt.gz')
        url = f'{STATIC_ASSET}{getattr(self,table)}/{taxon}.{getattr(self,table)}.txt.gz'
        if not file.is_file():
            print('Downloading...')
            with open(file, 'wb') as f:
                with httpx.stream("GET", url) as res:
                    for chunk in res.iter_raw():
                        f.write(chunk)
        return file


class Meta:
    def __init__(self, taxon: int = 9606, *, idents: Iterable):
        self.species = taxon
        self.symbols = idents
        self.cache = cache_path if cache_path else default_cache

    def map_id(self) -> pd.DataFrame:
        '''
        returned fields:
        ------
        queryIndex	position of the protein in your input (starting from position 0)
        stringId	STRING identifier
        ncbiTaxonId	NCBI taxon identifier
        taxonName	species name
        preferredName	common protein name
        annotation	protein annotation
        '''
        symbol = '\r'.join(self.symbols)
        res = client.post('/api/tsv/get_string_ids',
                data={'identifiers': symbol, 'species': self.species, 'limit': 1})
        string_ids = pd.read_csv(StringIO(res.text), sep='\t')
        print(f'Conversion rate is: {len(string_ids) / len(self.symbols)} !')
        return string_ids

    def map_id_local(self) -> pd.Series:
        pr_id_table = DbFile(STRING_VER).download(self.species, 'alias')
        protein_aliases: pd.DataFrame = pd.read_csv(pr_id_table, sep='\t')
        ## warning : may loss genes/proteins
        _to_stringid = protein_aliases.drop_duplicates(subset=['alias','#string_protein_id'])
        string_ids = _to_stringid[_to_stringid['alias'].isin(self.symbols)]\
                .groupby(['alias'])['#string_protein_id'].first()
        print(f'Conversion rate is: {len(string_ids) / len(self.symbols)} !')
        return string_ids

    def protein_info_local(self):
        pr_anno_table = Path.joinpath(self.cache, f'{self.species}.{DbFile(STRING_VER).info}.txt.gz')
        if not pr_anno_table.is_file():
            DbFile(STRING_VER).download(self.species, 'info')
        protein_anno: pd.DataFrame = pd.read_csv(pr_anno_table, sep='\t')
        pr_anno = protein_anno[(protein_anno['preferred_name']==self.symbols) 
                | (protein_anno['#string_protein_id']==self.symbols)]
        #return pr_df[['#string_protein_id', 'protein_size', 'annotation']]
        return pr_anno
    

if __name__ == '__main__':
    def main():
        init()
        DbFile(STRING_VER).download(taxon= sys.argv[1], table=sys.argv[2])
