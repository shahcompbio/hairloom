# Hairloom

**Hairloom** is a Python package for analyzing split-read alignments from long-read sequencing data. It provides tools to extract read-level split-read tables from BAM files and process these alignments into breakpoint, segment, and translocation tables for downstream analysis.

## Features

- Extract and normalize split-read alignments from BAM files.
- Generate tables for:
    - **Breakpoints**
    - **Genomic Segments**
    - **Translocations**
- Designed for compatibility with long-read sequencing data from platforms like ONT and PacBio.
- Simple integration with downstream genomic workflows.

---

## Installation

```bash
pip install hairloom
```

---

## Usage

Below is an example workflow using **Hairloom**:

### Input

- A BAM file (`example.bam`).
- A specific genomic region (`chrom`, `start`, `end`).

### Step 1: Extract Read-Level Split Read Table

```python
import pysam
from hairloom import extract_read_data

# Open the BAM file
bam_file = pysam.AlignmentFile("example.bam", "rb")

# Extract split-read alignment table for a region
chrom, start, end = "chr1", 100000, 101000
read_table = extract_read_data(bam_file, chrom, start, end)

print(read_table.head())
```

**Output**: A table of split-read alignments with columns like:

- `qname`: Read name.
- `chrom`: Chromosome.
- `start`: Alignment start position.
- `end`: Alignment end position.
- `strand`: Strand information.
- `clip1`, `clip2`: Soft/hard clip lengths.
- `match`: Number of matched bases.
- `pclip1`: Strand-corrected clip length.

---

### Step 2: Generate Breakpoint Table

```python
from hairloom import make_brk_table, make_brk_supports

# Generate breakpoint supports
brk_supports = make_brk_supports([read_table])

# Create a breakpoint table
breakpoint_table = make_brk_table([read_table], brk_supports)

print(breakpoint_table.head())
```

**Output**: A table of breakpoints with columns:

- `chrom`: Chromosome.
- `pos`: Position.
- `ori`: Orientation (`+` or `-`).
- `support`: Support count.

---

### Step 3: Generate Segment Table

```python
from hairloom import make_seg_table

# Generate segment table
segment_supports = {(chrom, 100100, 100500): 10}  # Example segment supports
segment_table = make_seg_table([read_table], segment_supports)

print(segment_table.head())
```

**Output**: A table of genomic segments with columns:

- `chrom`: Chromosome.
- `pos1`: Start position.
- `pos2`: End position.
- `support`: Support count.

---

### Step 4: Generate Translocation Table

```python
from hairloom import make_tra_table

# Example translocation supports
tra_supports = {
    (("chr1", 100100, "+"), ("chr2", 200200, "-")): 5
}

# Create a translocation table
translocation_table = make_tra_table([read_table], tra_supports)

print(translocation_table.head())
```

**Output**: A table of translocations with columns:

- `chrom1`, `pos1`, `ori1`: First breakpoint information.
- `chrom2`, `pos2`, `ori2`: Second breakpoint information.
- `support`: Support count.

---

## Examples

### Full Workflow

```python
import pysam
from hairloom import (
    extract_read_data,
    make_brk_table,
    make_brk_supports,
    make_seg_table,
    make_tra_table,
)

# Open BAM file
bam_file = pysam.AlignmentFile("example.bam", "rb")

# Step 1: Extract read-level split-read table
chrom, start, end = "chr1", 100000, 101000
read_table = extract_read_data(bam_file, chrom, start, end)

# Step 2: Generate breakpoint table
brk_supports = make_brk_supports([read_table])
breakpoint_table = make_brk_table([read_table], brk_supports)

# Step 3: Generate segment table
segment_supports = {(chrom, 100100, 100500): 10}
segment_table = make_seg_table([read_table], segment_supports)

# Step 4: Generate translocation table
tra_supports = {
    (("chr1", 100100, "+"), ("chr2", 200200, "-")): 5
}
translocation_table = make_tra_table([read_table], tra_supports)

# Print results
print("Read Table:")
print(read_table.head())
print("\nBreakpoint Table:")
print(breakpoint_table.head())
print("\nSegment Table:")
print(segment_table.head())
print("\nTranslocation Table:")
print(translocation_table.head())
```

---

## Contributing

Contributions are welcome! Please submit issues and pull requests via the [GitHub repository](https://github.com/shahcompbio/hairloom).

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

Feel free to adapt this `README.md` to fit your exact requirements or additional features in the **Hairloom** package!
