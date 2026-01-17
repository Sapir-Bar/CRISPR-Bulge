import argparse
import pandas as pd
from Bio import Align


def ungap(x) -> str:
    return "" if pd.isna(x) else str(x).replace("-", "").strip()


def _restore_N_in_gapped_sg(gapped_sg: str, sg_original_ungapped: str) -> str:
    """Put 'N' back into the aligned sgRNA string at positions where original sgRNA had 'N'."""
    n_pos = {i for i, ch in enumerate(sg_original_ungapped) if ch.upper() == "N"}
    out = []
    i_ungapped = 0
    for ch in gapped_sg:
        if ch == "-":
            out.append("-")
        else:
            out.append("N" if i_ungapped in n_pos else ch)
            i_ungapped += 1
    return "".join(out)


def find_alignment(sg: str, ot: str, max_score: int = 7):
    """
    Returns: (aln_ot, aln_sg, mismatches, bulges, score) or None if score > max_score.
    mismatch weight = 1, bulge (gap char) weight = 3
    Supports 'N' in sgRNA by trying 4 substitutions (A/C/G/T) and picking best score.
    """
    sg0 = ungap(sg).upper()
    ot = ungap(ot).upper()
    if not sg0 or not ot:
        return None

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 0
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -3

    # NEW: expand sgRNA if it contains N (replace ALL Ns with same base)
    candidates = [sg0] if "N" not in sg0 else [sg0.replace("N", b) for b in "ACGT"]

    best = None  # (score, a_ot, a_sg, mm, bul)
    best_sg_aligned = None

    for sg_var in candidates:
        aln = aligner.align(sg_var, ot)[0]
        a_sg, a_ot = str(aln[0]), str(aln[1])

        mm = sum((x != y) and (x != "-") and (y != "-") for x, y in zip(a_sg, a_ot))
        bul = sum((x == "-") or (y == "-") for x, y in zip(a_sg, a_ot))
        score = mm + 3 * bul

        if (best is None) or (score < best[0]):
            best = (score, a_ot, a_sg, mm, bul)

    if best is None:
        return None

    score, a_ot, a_sg, mm, bul = best
    if score > max_score:
        return None

    # NEW: restore N in displayed aligned sgRNA
    a_sg = _restore_N_in_gapped_sg(a_sg, sg0)

    return (a_ot, a_sg, mm, bul, score)


def add_cols(df: pd.DataFrame, sg_col: str, ot_col: str, prefix: str, max_score: int):
    ot_aln, sg_aln, mm, bul, sc = [], [], [], [], []
    for sg, ot in zip(df[sg_col], df[ot_col]):
        r = find_alignment(sg, ot, max_score=max_score)
        if r is None:
            ot_aln.append(""); sg_aln.append(""); mm.append(""); bul.append(""); sc.append("")
        else:
            a_ot, a_sg, m, b, s = r
            ot_aln.append(a_ot); sg_aln.append(a_sg); mm.append(m); bul.append(b); sc.append(s)

    df[f"{prefix}Align.off-target"] = ot_aln
    df[f"{prefix}Align.sgRNA"] = sg_aln
    df[f"{prefix}Align.#Mismatches"] = mm
    df[f"{prefix}Align.#Bulges"] = bul
    df[f"{prefix}Align.Score"] = sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--max-score", type=int, default=20)
    ap.add_argument("--sg-col", default="on_target")
    ap.add_argument("--h1-col", default="h1")
    ap.add_argument("--h2-col", default="h2")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    for c in [args.sg_col, args.h1_col, args.h2_col]:
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")

    add_cols(df, args.sg_col, args.h1_col, "h1.", args.max_score)
    add_cols(df, args.sg_col, args.h2_col, "h2.", args.max_score)

    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()