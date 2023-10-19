from httpx import Client, HTTPTransport

client = Client(base_url="https://string-db.org/api",
                headers={'User-Agent': "Mozilla/5.0 (X11; Linux x86_64; rv:115.0) Gecko/20100101 Firefox/115.0"},
                follow_redirects=True, timeout=10,
                transport=HTTPTransport(retries=3))


string_versions = client.get("/json/version").json()
STRINGdb_latest = string_versions[0]
latest_api: str = STRINGdb_latest["stable_address"]
latest_ver: str = STRINGdb_latest["string_version"]


def choose_version(version: float=latest_ver) -> str:
    for ver in string_versions:
        if ver['string_version'] == str(version):
            return ver['stable_address']
    print(f"Your requested version is not supported now, try using {latest_ver} instead.")
    old_ver = str(version).replace('.','-')
    return f'https://version-{old_ver}.string-db.org'

STATIC_ASSET = 'https://stringdb-downloads.org/download/'
