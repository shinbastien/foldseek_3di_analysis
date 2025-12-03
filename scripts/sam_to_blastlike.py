#!/usr/bin/env python3
"""
Convert SAM alignments to BLAST-like alignment blocks using a reference FASTA.

Usage:
  python3 sam_to_blastlike.py reference.fasta in.sam > out.blast

If in.sam is '-' or omitted the script reads from stdin.

Notes:
- This script reconstructs the reference subsequence using RNAME and POS fields
  from the provided FASTA file and the CIGAR string. It supports CIGAR ops: =, X, M, I, D, S, H.
- The output is a minimal BLAST-like block similar to ssw_test's default (non -s) output.
"""
import sys, re, textwrap

def load_fasta(path):
    seqs = {}
    name = None
    buf = []
    fh = open(path, 'r', encoding='utf-8')
    for line in fh:
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(buf).replace('\n','')
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf).replace('\n','')
    fh.close()
    return seqs

cig_re = re.compile(r'(\d+)([MIDNSHP=X])')

def expand_alignment(ref_seq, pos1, cigar, query_seq):
    # pos1: 1-based leftmost reference coordinate
    rpos = pos1 - 1  # 0-based index into ref_seq
    qpos = 0
    aligned_ref = []
    aligned_q   = []
    consumed_ref = 0
    consumed_q = 0

    for length_s, op in cig_re.findall(cigar):
        length = int(length_s)
        if op in ('=', 'X', 'M'):
            for i in range(length):
                if rpos < 0 or rpos >= len(ref_seq) or qpos >= len(query_seq):
                    rbase = ref_seq[rpos] if 0 <= rpos < len(ref_seq) else 'N'
                    qbase = query_seq[qpos] if qpos < len(query_seq) else 'N'
                else:
                    rbase = ref_seq[rpos]
                    qbase = query_seq[qpos]
                aligned_ref.append(rbase)
                aligned_q.append(qbase)
                rpos += 1; qpos += 1
                consumed_ref += 1; consumed_q += 1
        elif op == 'I':
            # insertion to reference: consume query only
            for i in range(length):
                qbase = query_seq[qpos] if qpos < len(query_seq) else 'N'
                aligned_ref.append('-')
                aligned_q.append(qbase)
                qpos += 1
                consumed_q += 1
        elif op == 'D':
            # deletion from reference: consume reference only
            for i in range(length):
                rbase = ref_seq[rpos] if rpos < len(ref_seq) else 'N'
                aligned_ref.append(rbase)
                aligned_q.append('-')
                rpos += 1
                consumed_ref += 1
        elif op in ('S', 'H'):
            # soft/hard clip: consume query (S) or ignore (H)
            if op == 'S':
                qpos += length
                consumed_q += length
            # H: hard clip, no sequence present in SAM
        else:
            # unknown op: skip
            pass

    # compute coordinates: target begin = pos1, end = pos1 + consumed_ref - 1
    tbegin = pos1
    tend = pos1 + consumed_ref - 1 if consumed_ref > 0 else pos1

    # query begin/end: find first/last non-clipped aligned positions in query
    # We approximate query_begin as 1-based index of first consumed query base that was aligned (not S/H)
    # This approximation may differ from original ssw output if CIGAR is complex.
    q_aligned_indices = []
    q_index = 0
    for length_s, op in cig_re.findall(cigar):
        length = int(length_s)
        if op == 'S':
            q_index += length
        elif op in ('=', 'X', 'M'):
            for i in range(length):
                q_aligned_indices.append(q_index + 1)
                q_index += 1
        elif op == 'I':
            for i in range(length):
                q_aligned_indices.append(q_index + 1)
                q_index += 1
        elif op == 'D':
            # consumes reference only
            pass
        else:
            pass

    qbegin = q_aligned_indices[0] if q_aligned_indices else 1
    qend = q_aligned_indices[-1] if q_aligned_indices else 0

    return (''.join(aligned_ref), ''.join(aligned_q), tbegin, tend, qbegin, qend)

def make_midline(a_ref, a_q):
    mid = []
    for r, q in zip(a_ref, a_q):
        if r == q and r != '-' and q != '-':
            mid.append('|')
        elif r == '-' or q == '-':
            mid.append(' ')
        else:
            mid.append('*')
    return ''.join(mid)

def format_block(rname, qname, AS, ZS, tbegin, tend, qbegin, qend, a_ref, a_q, width=60):
    out = []
    out.append(f"target_name: {rname}")
    out.append(f"query_name: {qname}")
    out.append(f"optimal_alignment_score: {AS if AS is not None else 'NA'}\tsub-optimal_alignment_score: {ZS if ZS is not None else 'NA'}\tstrand: +\ttarget_begin: {tbegin}\ttarget_end: {tend}\tquery_begin: {qbegin}\tquery_end: {qend}")
    out.append('')

    # wrap and print blocks
    for i in range(0, len(a_ref), width):
        chunk_ref = a_ref[i:i+width]
        chunk_mid = make_midline(chunk_ref, a_q[i:i+width])
        chunk_q   = a_q[i:i+width]
        # compute local positions for pretty printing
        # For simplicity we won't compute exact end positions per wrapped chunk here; show original begin/end for whole
        out.append(f"Target:      {tbegin + i:>6}    {chunk_ref}    {tend if i+width>=len(a_ref) else ''}")
        out.append(f"             { ' ' * 6 }    {chunk_mid}")
        out.append(f"Query:        {qbegin + i:>6}    {chunk_q}    {qend if i+width>=len(a_q) else ''}")
        out.append('')
    return '\n'.join(out)

def parse_opt_fields(fields):
    AS = None; ZS = None
    for f in fields[11:]:
        if f.startswith('AS:i:'):
            try: AS = int(f.split(':')[-1])
            except: pass
        elif f.startswith('ZS:i:'):
            try: ZS = int(f.split(':')[-1])
            except: pass
    return AS, ZS

def main():
    if len(sys.argv) < 2:
        print('Usage: sam_to_blastlike.py reference.fasta [in.sam|-]', file=sys.stderr)
        sys.exit(2)
    ref_fasta = sys.argv[1]
    infpath = sys.argv[2] if len(sys.argv) > 2 else '-'

    refs = load_fasta(ref_fasta)

    if infpath == '-' or infpath == None:
        fh = sys.stdin
    else:
        fh = open(infpath, 'r', encoding='utf-8', errors='ignore')

    for line in fh:
        if line.startswith('@'): continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11: continue
        qname = fields[0]
        flag = int(fields[1]) if fields[1].isdigit() else 0
        rname = fields[2]
        pos = int(fields[3])
        cigar = fields[5]
        seq = fields[9]
        AS, ZS = parse_opt_fields(fields)

        if rname == '*' or cigar == '*':
            # unmapped
            continue

        if rname not in refs:
            print(f"# WARNING: reference {rname} not found in {ref_fasta}", file=sys.stderr)
            continue

        ref_seq = refs[rname]
        a_ref, a_q, tbegin, tend, qbegin, qend = expand_alignment(ref_seq, pos, cigar, seq)
        block = format_block(rname, qname, AS, ZS, tbegin, tend, qbegin, qend, a_ref, a_q)
        print(block)

    if fh is not sys.stdin:
        fh.close()

if __name__ == '__main__':
    main()
