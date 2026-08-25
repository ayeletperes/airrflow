#!/usr/bin/env python3
# Written by Ayelet Peres and released under the MIT license.
"""Restrict a germline reference_base tree to a single locus.

Selecting e.g. IGH means the reference must contain only IGH genes: leaving
IGKV/IGLV in place would let IgBLAST assign them, which is the whole point of
asking for a single locus. So the restricted tree holds only the wanted chains,
and the BLAST databases are rebuilt from it (bin/ref2igblast.sh) so no gene
outside the locus survives in a database.

The wanted files are hardlinked, not copied, and the unwanted ones are never
read: a reference tree is ~2 MB and every (species, locus) key would otherwise
duplicate all of it, which multiplies once personal germline sets add a key per
subject.

D genes are the exception: IGK, IGL, TRA and TRG have none, and AssignGenes.py's
--ddb is not optional, so the whole class D set is always kept. An IGK-restricted
run therefore sees exactly the D genes an IG run already gives it.
"""
import argparse
import os
import re
import shutil
import sys
from pathlib import Path

CHAIN = re.compile(r"(IG[HKL]|TR[ABGD])([VDJC])", re.IGNORECASE)
CLASS_D = {"IG": {"IGHD"}, "TR": {"TRBD", "TRDD"}}


def keep_chains(locus):
    """Chains to retain when restricting to `locus`, e.g. IGH -> IGHV/J/C + IGHD."""
    locus = locus.upper()
    if locus[:2] not in CLASS_D:
        sys.exit(f"Unsupported locus '{locus}': expected IGH, IGK, IGL, TRA, TRB, TRG or TRD.")
    return {f"{locus}V", f"{locus}J", f"{locus}C"} | CLASS_D[locus[:2]]


def chain_of(path):
    """The chain a reference FASTA holds, from its name, or None if unrecognised."""
    m = CHAIN.search(path.stem)
    return (m.group(1) + m.group(2)).upper() if m else None


def restrict(src, dest, locus, species=None):
    """Link `src` into `dest`, keeping only `species` and chains within `locus`."""
    keep = keep_chains(locus)
    src = Path(src)
    dropped = []

    def ignore(current, names):
        skip = set()
        # The reference ships every species; a reference keyed to one of them has
        # no use for the others, and building their databases is pure waste.
        if species and Path(current).resolve() == src.resolve():
            skip |= {n for n in names
                     if (src / n).is_dir() and n.lower() != species.lower()}
        for name in names:
            chain = chain_of(Path(name))
            # None means it is not a per-chain file (IMGT.yaml, a directory);
            # those are carried over untouched.
            if chain is not None and chain not in keep:
                skip.add(name)
                dropped.append(chain)
        return skip

    try:
        shutil.copytree(src, dest, ignore=ignore, copy_function=os.link)
    except OSError:
        # Hardlinks cannot cross filesystems; fall back to copying.
        shutil.rmtree(dest, ignore_errors=True)
        dropped.clear()
        shutil.copytree(src, dest, ignore=ignore)

    kept = [c for c in (chain_of(p) for p in sorted(dest.rglob("*.fasta"))) if c]
    return kept, dropped


def self_test():
    import tempfile
    assert keep_chains("IGH") == {"IGHV", "IGHJ", "IGHC", "IGHD"}
    assert keep_chains("IGK") == {"IGKV", "IGKJ", "IGKC", "IGHD"}, "IGK has no D of its own"
    assert keep_chains("TRB") == {"TRBV", "TRBJ", "TRBC", "TRBD", "TRDD"}
    assert chain_of(Path("imgt_human_IGHV.fasta")) == "IGHV"
    assert chain_of(Path("imgt_aa_human_IGKV.fasta")) == "IGKV"
    assert chain_of(Path("IMGT.yaml")) is None
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "ref" / "human" / "vdj"
        src.mkdir(parents=True)
        for chain in ("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV"):
            (src / f"imgt_human_{chain}.fasta").write_text(">x\nACGT\n")
        mouse = Path(tmp) / "ref" / "mouse" / "vdj"
        mouse.mkdir(parents=True)
        (mouse / "imgt_mouse_IGHV.fasta").write_text(">m\nACGT\n")
        (src.parent.parent / "IMGT.yaml").write_text("v: 1\n")
        out = Path(tmp) / "out"
        kept, dropped = restrict(Path(tmp) / "ref", out, "IGH", "human")
        assert not (out / "mouse").exists(), "other species must not be carried over"
        assert sorted(kept) == ["IGHD", "IGHJ", "IGHV"], kept
        assert sorted(dropped) == ["IGKJ", "IGKV", "IGLV"], dropped
        assert (out / "IMGT.yaml").exists(), "non-chain files must survive"
        assert not (out / "human" / "vdj" / "imgt_human_IGKV.fasta").exists()
        kept_file = out / "human" / "vdj" / "imgt_human_IGHV.fasta"
        assert kept_file.stat().st_nlink > 1, "kept files must be hardlinked, not copied"
    print("self-test ok")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-r", "--reference", type=Path, help="reference_base to restrict")
    p.add_argument("-l", "--locus", help="IGH, IGK, IGL, TRA, TRB, TRG or TRD")
    p.add_argument("-s", "--species", help="keep only this species' tree")
    p.add_argument("-o", "--out", type=Path, default=Path("reference_base"))
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()

    if a.self_test:
        return self_test()
    if not a.reference or not a.locus:
        p.error("-r and -l are required")

    kept, dropped = restrict(a.reference, a.out, a.locus, a.species)
    if not kept:
        sys.exit(f"No {a.locus} chain FASTAs found under {a.reference}; nothing to restrict to.")
    print(f"{a.locus}: kept {sorted(set(kept))}, dropped {sorted(set(dropped))}")


if __name__ == "__main__":
    main()
