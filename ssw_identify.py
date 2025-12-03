#!/usr/bin/env python3
# Usage:
#   python3 ssw_identity.py < in.sam > out.id.tsv
#   python3 ssw_identity.py in.sam > out.id.tsv
#
# 출력 컬럼:
# qname  rname  pos  cigar  L  NM  matches  identity_pct  AS  ZS
import sys, re

def parse_nm_and_score(fields):
    nm = None; AS = None; ZS = None
    for f in fields[11:]:
        if f.startswith("NM:i:"):
            try: nm = int(f.split(":")[-1])
            except: pass
        elif f.startswith("AS:i:"):
            try: AS = int(f.split(":")[-1])
            except: pass
        elif f.startswith("ZS:i:"):
            try: ZS = int(f.split(":")[-1])
            except: pass
    return nm, AS, ZS

cig_re = re.compile(r'(\d+)([MIDNSHP=XB])')

def compute_L_XID(cigar):
    # L = 합산(=,X,I,D,M)  (S/H/N/P는 정렬열에 미포함)
    # X,I,D 카운트도 반환(백업용)
    L = 0; Xc=Ic=Dc=Mc=Eqc=0
    for n,op in cig_re.findall(cigar):
        n = int(n)
        if op in ('=','X','I','D','M'):
            L += n
            if op == 'X': Xc += n
            elif op == 'I': Ic += n
            elif op == 'D': Dc += n
            elif op == 'M': Mc += n
            elif op == '=': Eqc += n
    return L, Xc, Ic, Dc, Mc, Eqc

def process_line(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6: return None
    if fields[0].startswith('@'): return None
    qname = fields[0]
    rname = fields[2]
    pos   = fields[3]
    cigar = fields[5]

    nm, AS, ZS = parse_nm_and_score(fields)
    L, Xc, Ic, Dc, Mc, Eqc = compute_L_XID(cigar)

    # NM이 있으면 matches = L - NM (표준 SAM 규칙: NM=mismatches+I+D)
    # 없으면 (=,X 표기 있는 경우) nm≈X+I+D 로 대체, M만 쓰인 경우엔 신뢰도 낮음.
    if nm is None:
        nm = Xc + Ic + Dc  # 대체 추정
    matches = max(L - nm, 0)
    identity = (matches / L * 100.0) if L > 0 else 0.0

    # 탭 구분 TSV 출력
    return "\t".join([
        qname, rname, pos, cigar,
        str(L), str(nm), str(matches),
        f"{identity:.2f}",
        str(AS) if AS is not None else "NA",
        str(ZS) if ZS is not None else "NA"
    ])

def main():
    wrote_header = False
    def write_header():
        print("\t".join(["qname","rname","pos","cigar","L","NM","matches","identity_pct","AS","ZS"]))
    if len(sys.argv) > 1 and sys.argv[1] != '-':
        fh = open(sys.argv[1], "r", encoding="utf-8", errors="ignore")
    else:
        fh = sys.stdin
    for line in fh:
        if line.startswith("@"):
            continue
        res = process_line(line)
        if res:
            if not wrote_header:
                write_header(); wrote_header = True
            print(res)
    if fh is not sys.stdin:
        fh.close()

if __name__ == "__main__":
    main()
