#!/usr/bin/env python3
# Written by Ayelet Peres and released under the MIT license.
"""Build a germline reference for one sample group.

Restricts a reference_base tree to a single locus, and/or overlays a subject's
personal germline set onto it. Wanted files are hardlinked, unwanted ones never
read. Class D genes are always kept: IGK, IGL, TRA and TRG have none of their
own and AssignGenes.py's --ddb is not optional.
"""
import argparse
import os
import re
import shutil
import sys
from pathlib import Path

CHAIN = re.compile(r"(IG[HKL]|TR[ABGD])([VDJCL])", re.IGNORECASE)
CLASS_D = {"IG": {"IGHD"}, "TR": {"TRBD", "TRDD"}}
# Locus is the GGS directory name; the segment is the file name.
GGS_SEGMENTS = {"V_gapped_asc.fasta": "V", "D_asc.fasta": "D", "J_asc.fasta": "J"}


def keep_chains(locus):
    # D is the one cross-locus exception: IGK/IGL/TRA/TRG have none of their own
    # and AssignGenes.py --ddb is not optional.
    locus = locus.upper()
    if locus[:2] not in CLASS_D:
        sys.exit(f"Unsupported locus '{locus}': expected IGH, IGK, IGL, TRA, TRB, TRG or TRD.")
    return {f"{locus}V", f"{locus}J", f"{locus}C", f"{locus}L"} | CLASS_D[locus[:2]]


def chain_of(path):
    m = CHAIN.search(path.stem)
    return (m.group(1) + m.group(2)).upper() if m else None


def restrict(src, dest, locus=None, species=None):
    # Hardlinks, never copies: a reference tree is ~2 MB and every key would
    # otherwise duplicate all of it.
    keep = keep_chains(locus) if locus else None
    src = Path(src)
    dropped = []

    def ignore(current, names):
        skip = set()
        # The reference ships every species; a reference keyed to one of them has
        # no use for the others, and building their databases is pure waste.
        if species and Path(current).resolve() == src.resolve():
            skip |= {n for n in names
                     if (src / n).is_dir() and n.lower() != species.lower()}
        if keep is None:
            return skip
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


def ggs_chains(ggs):
    found = {}
    for sub in sorted(Path(ggs).iterdir()):
        if not sub.is_dir():
            continue
        for name, segment in GGS_SEGMENTS.items():
            if (sub / name).is_file():
                found[sub.name.upper() + segment] = sub / name
    return found


def missing_locus(chains, required):
    have = {c[:-1] for c in chains}
    return next((locus for locus in sorted(required) if locus not in have), None)


def apply_ggs(ggs, dest, species, keep=None):
    vdj = Path(dest) / species / "vdj"
    if not vdj.is_dir():
        sys.exit(f"No {species} tree at {vdj}; the reference has nothing to overlay a germline set onto.")

    replaced = []
    for chain, source in sorted(ggs_chains(ggs).items()):
        if keep is not None and chain not in keep:
            continue
        # Every prefix, not just imgt_: an airrc_ file left behind would put
        # generic alleles straight back into the database ref2igblast.sh builds.
        for stale in vdj.glob("*.fasta"):
            if chain_of(stale) == chain:
                stale.unlink()
        target = vdj / f"imgt_{species}_{chain}.fasta"
        try:
            os.link(source, target)
        except OSError:
            shutil.copy(source, target)
        replaced.append(chain)
    return replaced


def self_test():
    import tempfile
    assert keep_chains("IGH") == {"IGHV", "IGHJ", "IGHC", "IGHL", "IGHD"}
    assert keep_chains("IGK") == {"IGKV", "IGKJ", "IGKC", "IGKL", "IGHD"}, "IGK borrows class D"
    assert chain_of(Path("imgt_human_IGKL.fasta")) == "IGKL"
    assert chain_of(Path("IMGT.yaml")) is None

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        vdj = tmp / "ref" / "human" / "vdj"
        vdj.mkdir(parents=True)
        for chain in ("IGHV", "IGHD", "IGHJ", "IGKV", "IGLV"):
            (vdj / f"imgt_human_{chain}.fasta").write_text(">generic\nACGT\n")
        (tmp / "ref" / "mouse").mkdir()
        (tmp / "ref" / "IMGT.yaml").write_text("v: 1\n")

        out = tmp / "out"
        kept, dropped = restrict(tmp / "ref", out, "IGH", "human")
        assert sorted(kept) == ["IGHD", "IGHJ", "IGHV"], kept
        assert sorted(dropped) == ["IGKV", "IGLV"], dropped
        assert not (out / "mouse").exists(), "other species are not carried over"
        assert (out / "IMGT.yaml").exists(), "non-chain files survive"
        assert (out / "human" / "vdj" / "imgt_human_IGHV.fasta").stat().st_nlink > 1, "must hardlink"

        ggs = tmp / "ggs" / "IGH"
        ggs.mkdir(parents=True)
        (ggs / "V_gapped_asc.fasta").write_text(">personal\nAC.GT\n")
        assert apply_ggs(tmp / "ggs", out, "human") == ["IGHV"]
        text = (out / "human" / "vdj" / "imgt_human_IGHV.fasta").read_text()
        assert "personal" in text and "." in text, "GGS replaces the chain, gaps intact"
    print("self-test ok")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-r", "--reference", type=Path, help="reference_base to restrict")
    p.add_argument("-l", "--locus", help="IGH, IGK, IGL, TRA, TRB, TRG or TRD")
    p.add_argument("-s", "--species", help="keep only this species' tree")
    p.add_argument("-g", "--ggs", type=Path, help="personal germline set to overlay")
    p.add_argument("--require-loci", default="", help="comma-separated loci the GGS must provide")
    p.add_argument("--subject", default="", help="subject id, for the --require-loci error message")
    p.add_argument("-o", "--out", type=Path, default=Path("reference_base"))
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()

    if a.self_test:
        return self_test()
    if not a.reference or not (a.locus or a.ggs):
        p.error("-r and at least one of -l/-g are required")
    if a.ggs and not a.species:
        p.error("-s is required with -g")

    kept, dropped = restrict(a.reference, a.out, a.locus, a.species)
    if a.locus and not kept:
        sys.exit(f"No {a.locus} chain FASTAs found under {a.reference}; nothing to restrict to.")
    print(f"{a.locus or 'all loci'}: kept {sorted(set(kept))}, dropped {sorted(set(dropped))}")

    if a.ggs:
        root = Path(a.ggs)
        chains = ggs_chains(root)
        if not chains:
            sys.exit(f"--ggs_input: ggs_path for subject '{a.subject}' has no "
                     f"<locus>/*.fasta germline set: {a.ggs}")
        # For a zipped germline set this is the first chance to look inside, so
        # the coverage check the samplesheet does for a directory runs here too.
        gap = missing_locus(chains, [x for x in a.require_loci.split(",") if x])
        if gap:
            sys.exit(f"--ggs_input: subject '{a.subject}' has samples needing locus {gap}, "
                     f"which the germline set does not provide: {a.ggs}")
        replaced = apply_ggs(root, a.out, a.species, keep_chains(a.locus) if a.locus else None)
        print(f"ggs {a.subject}: replaced {sorted(replaced)}")


if __name__ == "__main__":
    main()
