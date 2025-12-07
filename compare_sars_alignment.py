from Bio import SeqIO

alignment_file = "../result/SARS_Aligned.fasta"

records = list(SeqIO.parse(alignment_file, "fasta"))

if len(records) != 2:
    raise ValueError(f"Expected 2 sequences, found {len(records)}")

seq1 = str(records[0].seq)
seq2 = str(records[1].seq)

if len(seq1) != len(seq2):
    raise ValueError("Aligned sequences are not the same length")

total_len = len(seq1)
diff_count = 0
diff_positions = []

for i, (b1, b2) in enumerate(zip(seq1, seq2), start=1):
    if b1 != b2:
        diff_count += 1
        if len(diff_positions) < 20:
            diff_positions.append((i, b1, b2))

print(f"Total alignment length: {total_len}")
print(f"Number of differing positions: {diff_count}")
print(f"Percentage difference: {diff_count / total_len * 100:.3f}%")

print("\nFirst few differences (Position, SARS-CoV-2 base, SARS-CoV Tor2 base):")
for pos, b1, b2 in diff_positions:
    print(pos, b1, b2)

