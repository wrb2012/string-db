from httpx import Client

client = Client(base_url="https://string-db.org/api",
                headers={'User-Agent': "Mozilla/5.0 (X11; Linux x86_64; rv:115.0) Gecko/20100101 Firefox/115.0"},
                timeout=10,
                follow_redirects=True)


STRINGdb_latest: dict = client.get("/json/version").json()[0]
api_domain: str = STRINGdb_latest["stable_address"]
STRING_VER: str = STRINGdb_latest["string_version"]

STATIC_FILES = 'https://stringdb-downloads.org/download/'
