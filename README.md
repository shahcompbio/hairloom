# Hairloom

<!-- index.rst content start -->

[![coverage](https://github.com/shahcompbio/hairloom/actions/workflows/check-code-coverage.yml/badge.svg?event=push)](https://github.com/shahcompbio/hairloom/actions/workflows/check-code-coverage.yml)
[![pytest CI](https://github.com/shahcompbio/hairloom/actions/workflows/check-install.yml/badge.svg)](https://github.com/shahcompbio/hairloom/actions/workflows/check-install.yml)

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

- Requires `python` >= 3.7 and `pysam` for smooth installation.
- Requires `pandas` >=1.1.0 to run CLI commands (as per `requirements.txt`).

```bash
# conda install -y -c bioconda r pysam # if you don't have pysam installed
pip install hairloom
```

Check integrity of the installation with `pytest`:
```bash
$ pytest
============================ test session starts =============================
platform linux -- Python 3.11.8, pytest-8.1.1, pluggy-1.4.0
rootdir: /data1/shahs3/users/chois7/projects/hairloom
configfile: pyproject.toml
plugins: cov-5.0.0, anyio-4.3.0
collected 29 items

tests/test_cli.py ....                                                 [ 13%]
tests/test_collect.py ........                                         [ 41%]
tests/test_datatypes.py ...........                                    [ 79%]
tests/test_utils.py ......                                             [100%]

============================= 29 passed in 3.43s =============================
```

---

## Usage

Below is an example usage of the **Hairloom** CLI:

### Input

- A BAM file (`tests/data/test_reads.bam`)..
- A specific genomic region (`chrom`, `start`, `end`).


### 1. Extract Split-Read Alignments
Use the `extract` command to extract split-read alignments from a BAM file for a specified region. 
The output is a table listing the aligned fragments, ordered according to their position within each read.

**Command:**

```bash
hairloom extract <bam> <chrom> <start> <end>
```

**Example:**

```bash
hairloom extract tests/data/test_reads.bam chr1 50 150
```

**Output:**

A TSV-format of split-read alignment data given to STDOUT:

```bash
qname   chrom   start   end     strand  clip1   match   clip2   pclip1
read1   chr1    101     201     +       0       100     300     0
read1   chr2    300     400     +       100     100     200     100
read1   chr2    700     800     -       100     100     200     200
read1   chr2    900     1000    +       300     100     0       300
```

### 2. Generate Breakpoint Table
Use the `breakpoints` command to generate a table of breakpoints from a BAM file.

**Command:**

```bash
hairloom breakpoints <bam> <chrom> <start> <end>
```

**Example:**

```bash
hairloom breakpoints tests/data/test_reads.bam chr1 50 150
```

**Output:**

A TSV-format STDOUT stream of breakpoints:

```bash
chrom   pos     ori     support
chr1    201     +       1
chr2    300     -       1
chr2    400     +       1
chr2    700     -       1
chr2    800     +       1
chr2    900     -       1
```


### 3. Generate Segment Table
Use the `segments` command to generate a table of breakpoints from a BAM file.

**Command:**

```bash
hairloom segments <bam> <chrom> <start> <end>
```

**Example:**

```bash
hairloom breakpoints tests/data/test_reads.bam chr1 50 150
```

**Output:**

A TSV-format STDOUT stream of genomic segments:

```bash
chrom   pos1    pos2    support
chr2    300     400     1
chr2    700     800     1
```



### 4. Generate SV Table
Use the `svs` command to generate a table of breakpoints from a BAM file.

**Command:**

```bash
hairloom svs <bam> <chrom> <start> <end>
```

**Example:**

```bash
hairloom svs tests/data/test_reads.bam chr1 50 150
```

**Output:**

A TSV-format STDOUT stream of breakpoint pairs:

```bash
chrom1  pos1    ori1    chrom2  pos2    ori2    support
chr1    201     +       chr2    300     -       1
chr2    400     +       chr2    800     +       1
chr2    700     -       chr2    900     -       1
```


## Python Workflow

Below is an example workflow using **Hairloom**:

### Input

- A BAM file (`tests/data/test_reads.bam`)..
- A specific genomic region (`chrom`, `start`, `end`).

### 1. Extract Read-Level Split Read Table

```python
import pysam
from hairloom import extract_read_data

# Open the BAM file
bam_path = "tests/data/test_reads.bam"
bam_file = pysam.AlignmentFile(bam_path, "rb")

# Extract split-read alignment table for a region
chrom, start, end = "chr1", 50, 150
read_table = extract_read_data(bam_file, chrom, start, end)

print("Extracted Read-Level Table:")
print(read_table)
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

```
Extracted Read-Level Split-Read Table:
   qname chrom  start   end strand  clip1  match  clip2  pclip1
0  read1  chr1    101   201      +      0    100    300       0
1  read1  chr2    300   400      +    100    100    200     100
2  read1  chr2    700   800      -    100    100    200     200
3  read1  chr2    900  1000      +    300    100      0     300
```
---

### 2. Generate Breakpoint Table

```python
from hairloom import make_bundle, make_brk_table, make_brk_supports

# Make bundled BreakpointChain from read table
bundle = make_bundle(read_table)

# Create a breakpoint table
breakpoint_table = make_brk_table(bundle)

print(breakpoint_table)
```

**Output**: A table of breakpoints with columns:

- `chrom`: Chromosome.
- `pos`: Position.
- `ori`: Orientation (`+` or `-`).
- `support`: Support count.

```
  chrom  pos ori  support
0  chr1  201   +        1
1  chr2  300   -        1
2  chr2  400   +        1
3  chr2  700   -        1
4  chr2  800   +        1
5  chr2  900   -        1
```

---

### 3. Generate Segment Table

```python
from hairloom import make_seg_table

# Generate segment table
segment_table = make_seg_table(bundle)

print(segment_table.head())
```

**Output**: A table of genomic segments with columns:

- `chrom`: Chromosome.
- `pos1`: Start position.
- `pos2`: End position.
- `support`: Support count.

```
  chrom  pos1  pos2  support
0  chr2   300   400        1
1  chr2   700   800        1
```

---

### 4. Generate Translocation Table

```python
from hairloom import make_tra_table

# Create a translocation table
translocation_table = make_tra_table(bundle)

print(translocation_table.head())
```

**Output**: A table of translocations with columns:

- `chrom1`, `pos1`, `ori1`: First breakpoint information.
- `chrom2`, `pos2`, `ori2`: Second breakpoint information.
- `support`: Support count.

```
  chrom1  pos1 ori1 chrom2  pos2 ori2  support
0   chr1   201    +   chr2   300    -        1
1   chr2   400    +   chr2   800    +        1
2   chr2   700    -   chr2   900    -        1
```

---

## Snippets

### Full Workflow

```python
import pysam
from hairloom import (
    extract_read_data,
    make_brk_table,
    make_seg_table,
    make_tra_table,
)

# Open BAM file
bam_file = pysam.AlignmentFile("tests/data/test_reads.bam", "rb")

# Step 1: Extract read-level split-read table
chrom, start, end = "chr1", 50, 150
read_table = extract_read_data(bam_file, chrom, start, end)
bundle = make_bundle(read_table)

# Step 2: Generate breakpoint table
breakpoint_table = make_brk_table(bundle)

# Step 3: Generate segment table
segment_table = make_seg_table(bundle)

# Step 4: Generate translocation table
translocation_table = make_tra_table(bundle)

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

### Normalize an SV Table

The `normalize_sv_table` function sorts and normalizes the breakpoint pairs in a structural variant (SV) table, ensuring consistent ordering for downstream analyses.

```python
import pandas as pd
from hairloom import normalize_sv_table

# Example SV table
sv_data = {
    "chromosome_1": ["chr1", "chr2", "chr1", "chr3"],
    "position_1": [200, 500, 300, 100],
    "strand_1": ["+", "-", "+", "-"],
    "chromosome_2": ["chr1", "chr2", "chr1", "chr2"],
    "position_2": [100, 400, 400, 200],
    "strand_2": ["-", "+", "-", "+"]
}
sv_table = pd.DataFrame(sv_data)

# Normalize the SV table
normalized_sv = normalize_sv_table(sv_table)

print("Original SV Table:")
print(sv_table)
print("\nNormalized SV Table:")
print(normalized_sv)
```

**Output:**

```
Original SV Table:
  chromosome_1  position_1 strand_1 chromosome_2  position_2 strand_2
0         chr1         200        +         chr1         100       -
1         chr2         500        -         chr2         400       +
2         chr1         300        +         chr1         400       -
3         chr3         100        -         chr2         200       +

Normalized SV Table:
  chromosome_1  position_1 strand_1 chromosome_2  position_2 strand_2
0         chr1         100        -         chr1         200       +
1         chr2         400        +         chr2         500       -
2         chr1         300        +         chr1         400       -
3         chr2         200        +         chr3         100       -
```

### Get Structural Variant Type
The `get_svtype` function determines the type of structural variant (SV) represented by a `BreakpointPair`. It supports common SV types: translocation (TRA), inversion (INV), duplication (DUP), and deletion (DEL).

```python
from hairloom import Breakpoint, BreakpointPair, get_svtype

# Example BreakpointPairs
tra1 = BreakpointPair(Breakpoint("chr1", 100, "+"), Breakpoint("chr2", 200, "-"))
inv = BreakpointPair(Breakpoint("chr1", 100, "+"), Breakpoint("chr1", 200, "+"))
del_sv = BreakpointPair(Breakpoint("chr1", 300, "+"), Breakpoint("chr1", 500, "-"))
dup = BreakpointPair(Breakpoint("chr1", 500, "-"), Breakpoint("chr1", 700, "+"))

# Get SV types
print("SV Types:")
print(f"TRA: {get_svtype(tra1)}")  # Output: TRA
print(f"INV: {get_svtype(inv)}")    # Output: INV
print(f"DEL: {get_svtype(del_sv)}")  # Output: DEL
print(f"DUP: {get_svtype(dup)}")    # Output: DUP
```

**Output:**

```
SV Types:
TRA: TRA
INV: INV
DEL: DEL
DUP: DUP
```
---

## Contributing

Contributions are welcome! Please submit issues and pull requests via the [GitHub repository](https://github.com/shahcompbio/hairloom).

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

Feel free to adapt this `README.md` to fit your exact requirements or additional features in the **Hairloom** package!
