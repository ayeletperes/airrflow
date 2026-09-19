#!/usr/bin/env python3
"""
Generate the tiny synthetic paired-end dataset behind the `test_rm` profile.

Mimics the R24 macaque layout: R1 headers carry the UMI and TRIM annotations
appended after the Illumina comment, e.g.
  @M00001:1:000000000-AAAAA:1:1101:10:1 1:N:0:1|UMI=ACGTAC|TRIM=ACGTACGGATCC

Each sample contains:
  - overlapping pairs with high quality      -> pass assembly and quality filter
  - overlapping pairs with low quality       -> pass assembly, fail FilterSeq
  - non-overlapping random pairs             -> fail assembly

Usage: python3 make_test_data.py [outdir]   (default: this directory)
"""
import gzip
import os
import random
import sys

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_pair(f1, f2, name, r1, r2, q1, q2, umi, trim):
    f1.write(f"@{name} 1:N:0:1|UMI={umi}|TRIM={trim}\n{r1}\n+\n{q1}\n")
    f2.write(f"@{name} 2:N:0:1\n{r2}\n+\n{q2}\n")


def make_sample(outdir, sample, seed, n_good, n_lowq, n_fail):
    rng = random.Random(seed)
    templates = [rand_seq(rng, 420) for _ in range(5)]
    read_len = 250
    p1 = os.path.join(outdir, f"{sample}_R1.fastq.gz")
    p2 = os.path.join(outdir, f"{sample}_R2.fastq.gz")
    with gzip.open(p1, "wt") as f1, gzip.open(p2, "wt") as f2:
        i = 0
        for kind, n in (("good", n_good), ("lowq", n_lowq), ("fail", n_fail)):
            for _ in range(n):
                i += 1
                name = f"M00001:1:000000000-AAAAA:1:1101:{i}:{seed}"
                umi = rand_seq(rng, rng.randint(4, 8))
                trim = umi + "ACGGATCC" + "GTACGATCGA"
                if kind == "fail":
                    r1, r2 = rand_seq(rng, read_len), rand_seq(rng, read_len)
                else:
                    t = rng.choice(templates)
                    r1, r2 = t[:read_len], revcomp(t[-read_len:])
                q = "+" if kind == "lowq" else "I"  # Phred 10 vs 40
                write_pair(f1, f2, name, r1, r2, q * read_len, q * read_len, umi, trim)
    return p1, p2


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)
    for k, sample in enumerate(("S1", "S2", "S3")):
        p1, p2 = make_sample(outdir, sample, seed=k + 1, n_good=150, n_lowq=30, n_fail=20)
        # S3 is left uncompressed: the real data is a mix of gzipped and plain fastq
        if sample == "S3":
            for p in (p1, p2):
                with gzip.open(p, "rb") as src, open(p[:-3], "wb") as dst:
                    dst.write(src.read())
                os.remove(p)


if __name__ == "__main__":
    main()
