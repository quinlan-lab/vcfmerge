"""
Merge two VCFs with the same samples, or a subset of the same
samples, for example, a VCF of structural variants and one of
small (GATK-called) variants
This will output the variants in sorted order if the input is sorted.
For this reason, it is best to sort the smaller SV VCF so that it has
the same chromosome ordering as the GATK VCF.
"""

from __future__ import print_function

import sys
import gzip
import re
from collections import OrderedDict

xopen = lambda f: gzip.open(f) if f.endswith(".gz") else open(f)

def get_header(path):
    _patt=re.compile("(\w+)=<ID=([^,]+),?.*")
    header = dict(infos=OrderedDict(), formats=OrderedDict(),
                  contigs=OrderedDict(), filters=OrderedDict(),
                  other=[], samples=[])

    for i, line in enumerate(xopen(path)):
        if i == 0: continue # skip VCF header
        if line[0] != "#": break
        line = line.rstrip('\r\n')
        if line.startswith("##INFO"):
            _, id = _patt.search(line).groups(0)
            header['infos'][id] = line
        elif line.startswith("##FORMAT"):
            _, id = _patt.search(line).groups(0)
            header['formats'][id] = line
        elif line.startswith("##FILTER"):
            _, id = _patt.search(line).groups(0)
            header['filters'][id] = line
        elif line.startswith("##contig"):
            _, id = _patt.search(line).groups(0)
            header['contigs'][id] = line
        elif line.startswith("#CHROM\t"):
            toks = line.split("\t")
            samples = toks[9:]
            header['samples'] = samples
        else:
            header['other'].append(line)
    return header

def combine_headers(h1, h2):
    h1_idxs = [(i, s) for i, s in enumerate(h1['samples']) if s in h2['samples']]
    h2_idxs = [(h2['samples'].index(s), s) for i, s in h1_idxs]
    h1['idxs'] = h1_idxs
    h2['idxs'] = h2_idxs
    header = dict(infos=OrderedDict(), formats=OrderedDict(),
                  contigs=OrderedDict(), filters=OrderedDict(),
                  other=[], samples=[])

    header['samples'] = [s[1] for s in h1_idxs]
    for k in ('formats', 'infos', 'contigs', 'filters'):
        header[k] = h1[k].copy()
        for id, info in h2[k].items():
            if id in header[k]:
                if header[k][id] != info:
                    print("WARNING: differing headers for %s" % id, file=sys.stderr)
                    if "=Float" in info and not "=Float" in header[k][id]:
                        print("    %s vs (using)%s" % (header[k][id], info), file=sys.stderr)
                        header[k][id] = info
                    else:
                        print("    %s vs (using)%s" % (info, header[k][id]), file=sys.stderr)
            else:
                header[k][id] = info
    header['other'] = list(h1['other'])
    for o in h2['other']:
        if o in header['other']: continue
        header['other'].append(o)
    return header

def header_str(header):
    lines = ["##fileformat=VCFv4.1"]
    lines.extend(header['other'])
    for k in ('formats', 'infos', 'contigs', 'filters'):
        lines.extend(header[k].values())
    hl = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    hl.extend(header['samples'])
    lines.append("\t".join(hl))
    return "\n".join(lines)

def all_unknown_or_ref(samples):
    ss = (s.split(":", 1)[0] for s in samples)
    return all(s in ('.', './.', '0/0', '.|.', '0|0') for s in ss)

def write_vcf_body(path, hdr, file_idx, remove_ref=False, file=sys.stdout):
    idxs = [h[0] for h in hdr['idxs']]
    skipped, tot = 0, 0
    for line in xopen(path):
        if line[0] == "#": continue
        toks = line.rstrip().split("\t")
        samples = toks[9:]
        tot += 1

        samples = [samples[i] for i in idxs]
        if remove_ref and all_unknown_or_ref(samples):
            skipped += 1
            continue
        toks = toks[:9] + samples
        if toks[3] == "N" and toks[4][0] == "<":
            toks[3] = "." # fix lumpy ref call of 'N'

        yield toks[0], int(toks[1]), file_idx, toks
    if remove_ref:
        print(">> skipped %d ref/unknown variants out of %d from %s"
                % (skipped, tot, path), file=sys.stderr)

def main(vcf1, vcf2, remove_ref=False):
    h1 = get_header(vcf1)
    h2 = get_header(vcf2)
    header = combine_headers(h1, h2)
    print(header_str(header))

    pq = []
    gen1 = write_vcf_body(vcf1, h1, 0, remove_ref)
    gen2 = write_vcf_body(vcf2, h2, 1, remove_ref)

    gens = [gen1, gen2]

    import heapq
    for g in gens:
        heapq.heappush(pq, next(g))

    last_chrom = None
    while pq:
        chrom, pos, i, toks = heapq.heappop(pq)

        # switch chroms, clear out the heap and start over.
        if last_chrom != chrom and last_chrom is not None:
            print("\t".join(toks))
            chrom, pos, i, toks = heapq.heappop(pq)
            print("\t".join(toks))
            for k in range(2):
                try:
                    heapq.heappush(pq, next(gens[k]))
                except StopIteration:
                    pass
            if pq == []:
                break
            chrom, pos, i, toks = heapq.heappop(pq)

        last_chrom = chrom

        try:
            heapq.heappush(pq, next(gens[i]))
        except StopIteration:
            pass
        except:
            raise
        print("\t".join(toks))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--remove-ref", action="store_true", default=False,
                   help="remove variants where all shared samples are either hom-ref or unknown")
    p.add_argument("vcf_a")
    p.add_argument("vcf_b")

    a = p.parse_args()
    main(a.vcf_a, a.vcf_b, a.remove_ref)
