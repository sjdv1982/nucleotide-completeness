"""Indicates which unique (nonredundant, filtered) conformers are exclusively in terminal fragments

Terminal fragments are those that contain any of the first two or last two nucleotides.
No matter if it is dinucleotides or trinucleotides, these are the first two and last two fragments
"""

import sys

lib = sys.argv[1]
assert lib in ("dinuc", "trinuc")
fraglen = 2 if lib == "dinuc" else 3
motif = sys.argv[2]

segment_file = f"nucleotide-fragments/segments.txt"
segments = {}
for l in open(segment_file).readlines():
    if not len(l.strip()):
        continue
    code, start, length = l.split()
    start, length = int(start), int(length)
    if code not in segments:
        segments[code] = []
    segments[code].append((start, length))

origin_file = f"nucleotide-fragments/{lib}/origin/{motif}.txt"

result = []
with open(origin_file) as f:
    for l in f.readlines():
        only_terminal = True
        for item in l.split("/"):
            if not item.strip():
                continue
            fields = item.split()
            assert len(fields) == 3
            code = fields[0]
            frag = int(fields[1])
            for segstart, seglength in segments[code]:
                if frag >= segstart and frag < segstart + seglength:
                    from_start = frag - segstart
                    from_end = segstart + seglength - frag - fraglen
                    if from_start in (0, 1) or from_end in (0, 1):
                        # is terminal
                        break
                    else:
                        only_terminal = False
                        break
            if not only_terminal:
                break
        result.append(only_terminal)

for l in result:
    print(int(l))
