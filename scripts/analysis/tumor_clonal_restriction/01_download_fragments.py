"""Download all 36 GSE302113 ATAC fragment files with size verification + resume.

Verifies each file against the server Content-Length; skips already-complete files;
retries transient failures. Writes a JSON status log the caller can read.
"""
import json, os, time, urllib.request, sys

RAW = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/raw/GSE302113"
os.makedirs(RAW, exist_ok=True)
jobs = json.load(open("frag_jobs.json"))
STATUS = "frag_download_status.json"

def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", 0))

def fetch(url, dest, expected, chunk=1 << 20):
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
        while True:
            b = r.read(chunk)
            if not b:
                break
            fh.write(b)
    got = os.path.getsize(tmp)
    if expected and got != expected:
        raise IOError(f"size mismatch {got} != {expected}")
    os.replace(tmp, dest)
    return got

status = {}
total = 0
for i, j in enumerate(jobs):
    dest = os.path.join(RAW, j["fname"])
    try:
        exp = remote_size(j["url"])
        if os.path.exists(dest) and os.path.getsize(dest) == exp:
            status[j["gsm"]] = {"ok": True, "bytes": exp, "skipped": True}
            total += exp
            continue
        for attempt in range(4):
            try:
                n = fetch(j["url"], dest, exp)
                status[j["gsm"]] = {"ok": True, "bytes": n}
                total += n
                break
            except Exception as e:
                if attempt == 3:
                    status[j["gsm"]] = {"ok": False, "error": str(e)}
                else:
                    time.sleep(5 * (attempt + 1))
    except Exception as e:
        status[j["gsm"]] = {"ok": False, "error": str(e)}
    json.dump({"done": i + 1, "of": len(jobs), "total_bytes": total, "status": status},
              open(STATUS, "w"), indent=1)
    print(f"[{i+1}/{len(jobs)}] {j['gsm']} {status[j['gsm']]}  cum={total/1e9:.1f} GB", flush=True)

ok = sum(1 for v in status.values() if v.get("ok"))
print(f"\nDONE: {ok}/{len(jobs)} ok, {total/1e9:.2f} GB total")
json.dump({"done": len(jobs), "of": len(jobs), "total_bytes": total, "status": status,
           "complete": True, "n_ok": ok}, open(STATUS, "w"), indent=1)
