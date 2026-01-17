#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <input.csv> <ref.fa> <output.csv> <padding>"
  echo "Example: $0 input.csv hg38.fa output.csv 25"
  exit 1
fi

INPUT="$1"
REF="$2"
OUT="$3"
PAD="${4:-25}"

# Ensure fasta index exists
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

python3 - <<'PY' "$INPUT" "$REF" "$OUT" "$PAD"
import sys, subprocess, tempfile
from pathlib import Path
import pandas as pd
from Bio import SeqIO

inp, ref, out, pad = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

df = pd.read_csv(inp)
req = ["#chr","start","end","strand","OT_name"]
miss = [c for c in req if c not in df.columns]
if miss:
    raise SystemExit(f"Missing required columns: {miss}")

# If your input start is 1-based, uncomment:
# df["start"] = df["start"].astype(int) - 1

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    bed = td/"regions.bed"
    genome = td/"ref.genome"
    winbed = td/"win.bed"
    fa = td/"win.fa"

    # genome file from .fai
    fai = Path(ref + ".fai")
    genome.write_text("".join(
        f"{line.split()[0]}\t{line.split()[1]}\n"
        for line in fai.read_text().splitlines() if line.strip()
    ))

    # CSV -> BED (name=OT_name for stable merge)
    bed_df = pd.DataFrame({
        "chr": df["#chr"].astype(str),
        "start": df["start"].astype(int),
        "end": df["end"].astype(int),
        "name": df["OT_name"].astype(str),
        "score": 0,
        "strand": df["strand"].astype(str),
    }).drop_duplicates()

    bed_df.to_csv(bed, sep="\t", header=False, index=False)

    # bedtools: sort -> slop -> sort -> getfasta
    s1 = subprocess.check_output(["bedtools","sort","-i",str(bed),"-g",str(genome)])
    p1 = subprocess.run(["bedtools","slop","-i","-","-g",str(genome),"-b",str(pad),"-s"],
                        input=s1, check=True, stdout=subprocess.PIPE)
    p2 = subprocess.run(["bedtools","sort","-i","-","-g",str(genome)],
                        input=p1.stdout, check=True, stdout=subprocess.PIPE)
    subprocess.run(["bedtools","getfasta","-fi",ref,"-bed","-","-fo",str(fa),"-name","-s"],
                   input=p2.stdout, check=True)

    # FASTA -> dict[OT_name]=sequence
    windows = {rec.id.strip(): str(rec.seq).upper() for rec in SeqIO.parse(str(fa), "fasta")}

df["SiteWindow"] = df["OT_name"].map(windows)

missing = df["SiteWindow"].isna().sum()
if missing:
    print(f"[WARN] {missing} rows missing SiteWindow (OT_name mismatch).", file=sys.stderr)

df.to_csv(out, index=False)
print(f"[DONE] Wrote: {out}")
PY
