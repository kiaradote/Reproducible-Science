"""Download the small public file that a few tests use (GEO GSM2230757, Baron pancreas donor 1) into testdata/."""
import hashlib
import urllib.request
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2230nnn/GSM2230757/suppl/GSM2230757_human1_umifm_counts.csv.gz"
dest = Path(__file__).resolve().parent / "testdata" / "GSM2230757_human1_umifm_counts.csv.gz"
dest.parent.mkdir(exist_ok=True)
if dest.exists():
    print("already present:", dest)
else:
    urllib.request.urlretrieve(URL, dest)
    print("downloaded:", dest)
print("sha256", hashlib.sha256(dest.read_bytes()).hexdigest())
