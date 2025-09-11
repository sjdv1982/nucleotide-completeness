"""
Analyzes "failures" (conformers > 1 A from any other conformer) in terms of:
- terminality
- worst PDBs
"""

import sys
import itertools
import numpy as np

lib = sys.argv[1]
assert lib in ("dinuc", "trinuc")
fraglen = 2 if lib == "dinuc" else 3

bases = ("A", "C")
motifs = ["".join(s) for s in itertools.product(bases, repeat=fraglen)]

all_origins = []
all_fits = []
all_terminality = []
allcodes = set()
for motif in motifs:
    fit_file = f"nucleotide-library/output/closest-fit/{lib}-{motif}.txt"
    fits = np.loadtxt(fit_file, dtype=float)[:, 1]
    terminality_file = f"output/{lib}-{motif}-terminality.txt"
    terminality = np.loadtxt(terminality_file, dtype=bool)
    origin_file = f"nucleotide-fragments/{lib}/origin/{motif}.txt"

    origins = []
    for l in open(origin_file).readlines():
        single_pdb = True
        single_code = None
        for item in l.split("/"):
            if not item.strip():
                continue
            fields = item.split()
            assert len(fields) == 3
            code = fields[0][:4]
            if single_code is None:
                single_code = code
            else:
                if code != single_code:
                    single_pdb = False
                    break
        if single_pdb and single_code:
            origins.append(single_code)
        else:
            origins.append(None)
        allcodes.add(code)

    assert len(fits) == len(terminality) == len(origins), motif

    all_fits.append(fits)
    all_terminality.append(terminality)
    all_origins.append(origins)

seqlen = {}
segment_file = f"nucleotide-fragments/segments.txt"
for l in open(segment_file).readlines():
    if not len(l.strip()):
        continue
    code, start, length = l.split()
    code = code[:4]
    start, length = int(start), int(length)
    if code not in seqlen:
        seqlen[code] = 0
    seqlen[code] += length


fits = np.concatenate(all_fits, dtype=float)
terminality = np.concatenate(all_terminality, dtype=bool)
origins = sum(all_origins, [])
assert len(fits) == len(terminality) == len(origins)

failure_mask = (fits > 1).astype(bool)
print(f"Number of conformers:                   {len(origins)}")
print("")
print(
    f"Number of >1A failures:                 {failure_mask.sum()} ({failure_mask.sum()/len(origins)*100:.1f} % of all)"
)
print(
    f"Number of terminal-only conformers:     {terminality.sum()} ({terminality.sum()/len(origins)*100:.1f} % of all)"
)
termfail = (terminality & failure_mask).sum()
print(f"Number of >1A terminal-only failures:   {termfail}")
print(f"    {termfail/terminality.sum()*100:.1f} % of terminal-only conformers")
print(f"    {termfail/failure_mask.sum()*100:.1f} % of failures")
print("")


def distr(origins, failure_mask, alltxt):
    print("Distribution of failures among most prevalent PDB codes")
    fail_origins = [ori for ori, fail in zip(origins, failure_mask) if fail]
    pdbs, counts = np.unique(fail_origins, return_counts=True)
    weighted_pdbcounts = {code: 0 for code in seqlen}
    pdbcounts = {code: 0 for code in seqlen}
    for k, v in zip(pdbs, counts):
        pdbcounts[k] = v
        weighted_pdbcounts[k] = v / seqlen[k]
    wpdbs = list(weighted_pdbcounts.keys())
    wcounts = [weighted_pdbcounts[code] for code in wpdbs]
    sorted_codes = [
        wpdbs[ind] for ind in np.argsort(wcounts)[::-1] if wpdbs[ind] is not None
    ]
    for pct in 1, 2, 5, 10:
        ntop = np.round(pct / 100 * len(allcodes)).astype(int)
        top_codes = set(sorted_codes[:ntop])
        pct_frag = len([ori for ori in origins if ori in top_codes])
        pct_fail = sum([pdbcounts[code] for code in top_codes])
        print(
            f"    Worst {pct} % of PDBs: {pct_fail}, {pct_fail/failure_mask.sum()*100:.1f} % of {alltxt} failures, {pct_frag/len(origins)*100:.1f} % of {alltxt} fragments"
        )
    notop_all = len([ori for ori in origins if ori not in top_codes])
    notop_fail = len(
        [n for n, ori in enumerate(origins) if failure_mask[n] and ori not in top_codes]
    )

    print(f"    Remaining 90 %: {notop_fail} failures among {notop_all} fragments")
    print(
        f"        {notop_fail/notop_all * 100:.1f} % of {alltxt} fragments in those structures"
    )
    print(f"        {notop_fail/failure_mask.sum() * 100:.1f} % of {alltxt} failures")
    print(
        f"        Those structures contain {notop_all/len(origins) * 100:.1f} % of {alltxt} fragments"
    )


distr(origins, failure_mask, "all")
print()
origins_nonterm = [ori for ori, term in zip(origins, terminality) if not term]
failure_mask_nonterm = failure_mask[~terminality]
distr(origins_nonterm, failure_mask_nonterm, "all non-terminal-only")
