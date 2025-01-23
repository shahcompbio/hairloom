import pysam

# Create reference sequences for chr1 and chr2
fasta_path = "test_ref.fasta"
with open(fasta_path, "w") as fasta_file:
    fasta_file.write(">chr1\n")
    fasta_file.write("A" * 1000 + "\n")  # 1000 As
    fasta_file.write(">chr2\n")
    fasta_file.write("T" * 1000 + "\n")  # 1000 Ts

# Index the FASTA file
pysam.faidx(fasta_path)

# Create the BAM file
bam_path = "test_reads.bam"
header = {
    "HD": {"VN": "1.0"},
    "SQ": [
        {"SN": "chr1", "LN": 1000},
        {"SN": "chr2", "LN": 1000},
    ],
}

# Open BAM file for writing
with pysam.AlignmentFile(bam_path, "wb", header=header) as bam_file:
    # Define 4 fragments for read1
    # Fragment 1: aligned to chr1
    read1_frag1 = pysam.AlignedSegment()
    read1_frag1.query_name = "read1"
    read1_frag1.query_sequence = "A" * 100
    read1_frag1.flag = 0
    read1_frag1.reference_id = 0  # chr1
    read1_frag1.reference_start = 100
    read1_frag1.mapping_quality = 60
    read1_frag1.cigar = [(0, 100), (5, 300)]  # 100M300H
    read1_frag1.tags = [("SA", "chr2,300,+,100H100M200H,60,100;chr2,700,-,100H100M200H,60,200;chr2,900,+,300H100M,60,300")]
    bam_file.write(read1_frag1)

    # Fragment 2: aligned to chr2
    read1_frag2 = pysam.AlignedSegment()
    read1_frag2.query_name = "read1"
    read1_frag2.query_sequence = "T" * 100
    read1_frag2.flag = 256
    read1_frag2.reference_id = 1  # chr2
    read1_frag2.reference_start = 300
    read1_frag2.mapping_quality = 60
    read1_frag2.cigar = [(5, 100), (0, 100), (5, 200)]  # 100H100M200H
    read1_frag2.tags = [("SA", "chr1,100,+,100M300H,60,0;chr2,700,-,100H100M200H,60,200;chr2,900,+,300H100M,60,300")]
    bam_file.write(read1_frag2)

    # Fragment 3: aligned to chr2
    read1_frag3 = pysam.AlignedSegment()
    read1_frag3.query_name = "read1"
    read1_frag3.query_sequence = "T" * 100
    read1_frag3.flag = 272
    read1_frag3.reference_id = 1  # chr2
    read1_frag3.reference_start = 700
    read1_frag3.mapping_quality = 60
    read1_frag3.cigar = [(5, 100), (0, 100), (5, 200)]  # 100H100M200H, but reverse mapped
    read1_frag3.tags = [("SA", "chr1,100,+,100M300H,60,0;chr2,300,+,100H100M200H,60,100;chr2,900,+,300H100M,60,300")]
    bam_file.write(read1_frag3)

    # Fragment 4: aligned to chr2
    read1_frag4 = pysam.AlignedSegment()
    read1_frag4.query_name = "read1"
    read1_frag4.query_sequence = "T" * 100
    read1_frag4.flag = 256
    read1_frag4.reference_id = 1  # chr2
    read1_frag4.reference_start = 900
    read1_frag4.mapping_quality = 60
    read1_frag4.cigar = [(5, 300), (0, 100)]  # 300H100M
    read1_frag4.tags = [("SA", "chr1,100,+,100M300H,60,0;chr2,300,+,100H100M200H,60,100;chr2,700,-,100M,60,200")]
    bam_file.write(read1_frag4)

# Index the BAM file
pysam.index(bam_path)

print(f"Test BAM file created: {bam_path}")
