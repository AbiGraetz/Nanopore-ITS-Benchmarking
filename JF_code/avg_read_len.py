#!/usr/bin/env python3
# avg_read_len.py

from pathlib import Path
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) != 2:
        print(f"usuage: {Path(sys.argv[0]).name} <input.fasta>", file=sys.stderr)
        sys.exit(1)

    fasta_path = Path(sys.argv[1])
    if not fasta_path.is_file():
        print(f"no file {fasta_path}", file=sys.stderr)
        sys.exit(1)

    total_len = 0
    seq_count = 0

    # 逐條讀取序列以節省記憶體
    for record in SeqIO.parse(fasta_path, "fasta"):
        total_len += len(record.seq)
        seq_count += 1

    if seq_count == 0:
        print("no read counts。")
        sys.exit(0)

    avg_len = total_len / seq_count
    print(f"read count：{seq_count}")
    print(f"avg length：{avg_len:.2f} bp")

if __name__ == "__main__":
    main()



