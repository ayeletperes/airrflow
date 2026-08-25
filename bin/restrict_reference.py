#!/usr/bin/env python3
# Written by Ayelet Peres and released under the MIT license.
"""Restrict a germline reference_base tree to a single locus, and/or overlay a
personal genomic germline set (GGS) onto it.

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
run therefore sees exactly the D genes an IG run already gives it. The same rule
holds for a GGS: only IGH ships a D file, so a GGS subject keeps the generic
class D set unless its own germline set replaces it.

A GGS (-g) is laid out as <ggs>/<LOCUS>/{V_gapped_asc,V_asc,D_asc,J_asc}.fasta --
the locus is the directory name, not part of the file name. The gapped V is the
one used: reference_base is IMGT-gapped and CreateGermlines depends on it;
ref2igblast.sh strips the gaps itself when it builds the BLAST databases. For
each chain the GGS provides, every generic FASTA for that chain is removed and
replaced, so no generic allele survives under any prefix.
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
    """Chains to retain when restricting to `locus`, e.g. IGH -> IGHV/J/C/L + IGHD.

    Leaders (L) follow their own locus: D is the one cross-locus exception, and
    it exists only because AssignGenes.py's --ddb is not optional.
    """
    locus = locus.upper()
    if locus[:2] not in CLASS_D:
        sys.exit(f"Unsupported locus '{locus}': expected IGH, IGK, IGL, TRA, TRB, TRG or TRD.")
    return {f"{locus}V", f"{locus}J", f"{locus}C", f"{locus}L"} | CLASS_D[locus[:2]]


def chain_of(path):
    """The chain a reference FASTA holds, from its name, or None if unrecognised."""
    m = CHAIN.search(path.stem)
    return (m.group(1) + m.group(2)).upper() if m else None


def restrict(src, dest, locus=None, species=None):
    """Link `src` into `dest`, keeping only `species` and chains within `locus`.

    `locus` of None keeps every chain: a GGS subject that asked for a whole
    receptor class restricts nothing, it only overlays its own alleles.
    """
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
    """{chain: fasta} for a GGS tree, e.g. IGH/V_gapped_asc.fasta -> IGHV."""
    found = {}
    for sub in sorted(Path(ggs).iterdir()):
        if not sub.is_dir():
            continue
        for name, segment in GGS_SEGMENTS.items():
            if (sub / name).is_file():
                found[sub.name.upper() + segment] = sub / name
    return found


def missing_locus(chains, required):
    """The first required locus the GGS has no chain for, or None."""
    have = {c[:-1] for c in chains}
    return next((locus for locus in sorted(required) if locus not in have), None)


def apply_ggs(ggs, dest, species, keep=None):
    """Replace the generic chains under `dest` with the ones the GGS provides."""
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
    assert keep_chains("IGK") == {"IGKV", "IGKJ", "IGKC", "IGKL", "IGHD"}, "IGK has no D of its own"
    assert keep_chains("TRB") == {"TRBV", "TRBJ", "TRBC", "TRBL", "TRBD", "TRDD"}
    assert chain_of(Path("imgt_human_IGHV.fasta")) == "IGHV"
    assert chain_of(Path("imgt_aa_human_IGKV.fasta")) == "IGKV"
    assert chain_of(Path("imgt_human_IGKL.fasta")) == "IGKL", "leader files are per-chain too"
    assert chain_of(Path("IMGT.yaml")) is None
    assert missing_locus({"IGHV", "IGHJ"}, ["IGH"]) is None
    assert missing_locus({"IGHV"}, ["IGH", "IGK", "IGL"]) == "IGK"

    def make_reference(root):
        vdj = root / "human" / "vdj"
        vdj.mkdir(parents=True)
        for chain in ("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV"):
            (vdj / f"imgt_human_{chain}.fasta").write_text(">generic\nACGT\n")
        mouse = root / "mouse" / "vdj"
        mouse.mkdir(parents=True)
        (mouse / "imgt_mouse_IGHV.fasta").write_text(">m\nACGT\n")
        leader = root / "human" / "leader"
        leader.mkdir(parents=True)
        for chain in ("IGHL", "IGKL"):
            (leader / f"imgt_human_{chain}.fasta").write_text(">generic\nACGT\n")
        (root / "IMGT.yaml").write_text("v: 1\n")
        return vdj

    def make_ggs(root, loci=("IGH", "IGK")):
        for locus in loci:
            (root / locus).mkdir(parents=True)
            for name in GGS_SEGMENTS:
                if name == "D_asc.fasta" and locus != "IGH":
                    continue
                (root / locus / name).write_text(f">{locus}-personal\nAC...GT\n")
        return root

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        make_reference(tmp / "ref")
        out = tmp / "out"
        kept, dropped = restrict(tmp / "ref", out, "IGH", "human")
        assert not (out / "mouse").exists(), "other species must not be carried over"
        assert sorted(kept) == ["IGHD", "IGHJ", "IGHL", "IGHV"], kept
        assert sorted(dropped) == ["IGKJ", "IGKL", "IGKV", "IGLV"], dropped
        assert (out / "human" / "leader" / "imgt_human_IGHL.fasta").exists()
        assert not (out / "human" / "leader" / "imgt_human_IGKL.fasta").exists(), \
            "a restricted reference must not keep another locus' leaders"
        assert (out / "IMGT.yaml").exists(), "non-chain files must survive"
        assert not (out / "human" / "vdj" / "imgt_human_IGKV.fasta").exists()
        kept_file = out / "human" / "vdj" / "imgt_human_IGHV.fasta"
        assert kept_file.stat().st_nlink > 1, "kept files must be hardlinked, not copied"

    with tempfile.TemporaryDirectory() as tmp:
        # No locus restriction: the GGS overlays IGH and IGK, IGL stays generic.
        tmp = Path(tmp)
        make_reference(tmp / "ref")
        (tmp / "ref" / "human" / "vdj" / "airrc_human_IGHV.fasta").write_text(">generic\nACGT\n")
        make_ggs(tmp / "ggs")
        out = tmp / "out"
        restrict(tmp / "ref", out, None, "human")
        replaced = apply_ggs(tmp / "ggs", out, "human")
        assert replaced == ["IGHD", "IGHJ", "IGHV", "IGKJ", "IGKV"], replaced
        vdj = out / "human" / "vdj"
        assert not (vdj / "airrc_human_IGHV.fasta").exists(), "every prefix of a replaced chain must go"
        assert "personal" in (vdj / "imgt_human_IGHV.fasta").read_text()
        assert "." in (vdj / "imgt_human_IGHV.fasta").read_text(), "IMGT gaps must survive"
        assert (vdj / "imgt_human_IGHV.fasta").stat().st_nlink > 1, "GGS files must be hardlinked"
        assert "generic" in (vdj / "imgt_human_IGLV.fasta").read_text(), "uncovered loci stay generic"

    with tempfile.TemporaryDirectory() as tmp:
        # IGK-restricted: the GGS IGKV/IGKJ are used and its IGHV/IGHJ must not
        # sneak IGH back in. IGHD is the class D set IGK borrows, so the
        # subject's own D replaces it, exactly as it would in an unrestricted run.
        tmp = Path(tmp)
        make_reference(tmp / "ref")
        make_ggs(tmp / "ggs")
        out = tmp / "out"
        restrict(tmp / "ref", out, "IGK", "human")
        replaced = apply_ggs(tmp / "ggs", out, "human", keep_chains("IGK"))
        assert replaced == ["IGHD", "IGKJ", "IGKV"], replaced
        vdj = out / "human" / "vdj"
        assert not (vdj / "imgt_human_IGHV.fasta").exists(), "restriction wins over the GGS"

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
