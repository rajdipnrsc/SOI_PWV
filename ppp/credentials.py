"""NASA Earthdata credentials for CDDIS downloads (TDS § 12.3: "Auth: .netrc for Earthdata (CDDIS)").

Looked up in this order (first one found wins):
  1. environment variables EARTHDATA_USERNAME / EARTHDATA_PASSWORD
  2. credentials.json in the working directory or next to pwv_ppp.py (git-ignored), written by
         python -m ppp.credentials
  3. settings.EARTHDATA_USERNAME / settings.EARTHDATA_PASSWORD (allowed, but never commit them)
  4. ~/.netrc entry "machine urs.earthdata.nasa.gov login USER password PASS" (used by requests automatically)

Credentials are sent only to the Earthdata hosts (settings.EARTHDATA_HOSTS), never written to logs, manifests or
the configuration hash.
"""
import getpass
import json
import os
import sys

from . import settings

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILE_NAME = "credentials.json"


def _candidates():
    return [os.path.abspath(FILE_NAME), os.path.join(HERE, FILE_NAME)]


def earthdata():
    """(username, password, source) or (None, None, source) if only ~/.netrc / nothing is available."""
    u, p = os.environ.get("EARTHDATA_USERNAME"), os.environ.get("EARTHDATA_PASSWORD")
    if u and p:
        return u, p, "environment"
    for f in _candidates():
        if os.path.exists(f):
            try:
                with open(f) as fh:
                    d = json.load(fh)
                if d.get("earthdata_username") and d.get("earthdata_password"):
                    return d["earthdata_username"], d["earthdata_password"], f
            except (OSError, ValueError):
                pass
    if settings.EARTHDATA_USERNAME and settings.EARTHDATA_PASSWORD:
        return settings.EARTHDATA_USERNAME, settings.EARTHDATA_PASSWORD, "settings.py"
    try:
        import netrc
        a = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
        if a:
            return None, None, "~/.netrc"
    except (OSError, netrc.NetrcParseError):
        pass
    return None, None, "none"


def main():
    """Interactive one-time setup: ask for the Earthdata login, test it on CDDIS, store it in credentials.json."""
    print("NASA Earthdata login (create one at https://urs.earthdata.nasa.gov/users/new)")
    u = input("Username: ").strip()
    p = getpass.getpass("Password: ")
    from .products import Downloader
    import tempfile
    dl = Downloader(tempfile.mkdtemp())
    dl.earthdata_override = (u, p)
    test = settings.EARTHDATA_TEST_URL
    ok = dl.get(test, os.path.join(dl.cache_dir, "test"))
    if not ok:
        print(f"Login test failed ({test}); see the messages above. Nothing was saved.")
        return 1
    path = os.path.join(HERE, FILE_NAME)
    with open(path, "w") as fh:
        json.dump({"earthdata_username": u, "earthdata_password": p}, fh)
    os.chmod(path, 0o600)
    print(f"Login works. Saved to {path} (readable only by you, ignored by git).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
