import pandas as pd

from .datatypes import Breakpoint, BreakpointChain

def enumerate_breakpoints(df):
    """
    Enumerates breakpoints from a DataFrame of genomic fragments.

    This function generates a `BreakpointChain` object containing `Breakpoint` objects
    based on the start and end positions of fragments in the input DataFrame. Breakpoints
    at the beginning and end of reads are treated differently to omit read boundaries.

    Args:
        df (pandas.DataFrame): A DataFrame containing fragment information, with the following columns:
            - 'chrom': The chromosome name (str).
            - 'start': The start position of the fragment (int).
            - 'end': The end position of the fragment (int).
            - 'strand': The strand information ('+' or '-').

    Returns:
        BreakpointChain: A chain of `Breakpoint` objects enumerated from the DataFrame.

    Notes:
        - For the first fragment in the DataFrame, only the end position is included as a breakpoint.
        - For the last fragment in the DataFrame, only the start position is included as a breakpoint.
        - For middle fragments, both start and end positions are included, with orientations determined
          by the strand.

    Example:
        >>> import pandas as pd
        >>> from your_module import Breakpoint, BreakpointChain, enumerate_breakpoints
        >>> df = pd.DataFrame({
        ...     'chrom': ['chr1', 'chr1', 'chr1'],
        ...     'start': [100, 200, 300],
        ...     'end': [150, 250, 350],
        ...     'strand': ['+', '+', '-']
        ... })
        >>> brks = enumerate_breakpoints(df)
        >>> print(brks)
        [chr1:150:+, chr1:200:-, chr1:250:+, chr1:300:-]
    """
    df = df.reset_index(drop=True)
    ix_start, ix_end = 0, df.shape[0] - 1
    brks = BreakpointChain([])
    for rix, row in df.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        if rix == ix_start: # omit read start from breakpoints
            ori = '+' if strand == '+' else '-'
            pos = end if strand == '+' else start
            brk = Breakpoint(chrom, pos, ori)
            _brks = [brk]
        elif rix == ix_end: # omit read end from breakpoints
            ori = '-' if strand == '+' else '+'
            pos = start if strand == '+' else end
            brk = Breakpoint(chrom, pos, ori)
            _brks = [brk]
        else:
            brk_start = Breakpoint(chrom, start, '-')
            brk_end = Breakpoint(chrom, end, '+')
            _brks = [brk_start, brk_end] if strand == '+' else [brk_end, brk_start]
        brks += _brks
    return brks

def get_secondaries(read):
    """
    Extracts secondary alignments from a sequencing read's 'SA' tag.

    This function retrieves secondary alignment information stored in the 'SA' tag
    of a sequencing read. It parses the tag and returns a list of secondary alignments.

    Args:
        read (pysam.AlignedSegment): A sequencing read object, typically from a 
            BAM or SAM file, that includes tags.

    Returns:
        list of str: A list of secondary alignment strings parsed from the 'SA' tag.
        If the 'SA' tag is not present, an empty list is returned.

    Example:
        >>> from pysam import AlignedSegment
        >>> read = AlignedSegment()
        >>> read.tags = [('SA', 'chr1,100,+,60M,60;chr2,200,-,60M,60;')]
        >>> secondaries = get_secondaries(read)
        >>> print(secondaries)
        ['chr1,100,+,60M,60', 'chr2,200,-,60M,60']

    Notes:
        - The 'SA' tag stores secondary alignments in a semicolon-separated format.
        - Each secondary alignment string typically contains the following fields:
            [chromosome, position, strand, CIGAR, mapping quality].
    """
    secondaries = [a[1].split(';') for a in read.tags if a[0] == 'SA']
    if len(secondaries) > 0:
        secondaries = secondaries[0]
        secondaries = list(filter(lambda x: x!='', secondaries))
        return secondaries
    return []

def make_split_read_table(alignments):
    """
    Creates a table summarizing split-read alignments.

    This function processes a list of alignment objects, extracting key attributes
    and organizing them into a structured pandas DataFrame. The resulting table
    is sorted and deduplicated based on read names and alignment positions.

    Args:
        alignments (list): A list of alignment objects. Each alignment object is expected
            to have the following attributes:
            - `read_name` (str): The query name of the read.
            - `refname` (str): The reference sequence name (chromosome or contig).
            - `start` (int): The start position of the alignment on the reference.
            - `end` (int): The end position of the alignment on the reference.
            - `strand` (str): The strand information ('+' or '-').
            - `clip1` (int): The length of the first soft/hard clip before the matched region.
            - `clip2` (int): The length of the second soft/hard clip after the matched region.
            - `match` (int): The total number of matched bases in the alignment.
            - `pclip1` (int): The strand-corrected length of the first clip.

    Returns:
        pandas.DataFrame: A DataFrame with the following columns:
            - `qname`: Query name of the read.
            - `chrom`: Reference sequence name (chromosome or contig).
            - `start`: Start position of the alignment on the reference.
            - `end`: End position of the alignment on the reference.
            - `strand`: Strand information ('+' or '-').
            - `clip1`: Length of the first soft/hard clip before the matched region.
            - `match`: Total number of matched bases in the alignment.
            - `clip2`: Length of the second soft/hard clip after the matched region.
            - `pclip1`: Strand-corrected length of the first clip.

    Notes:
        - The function removes duplicate rows based on the full set of columns.
        - The DataFrame is sorted by `qname` and `pclip1` for consistent organization.

    Example:
        >>> from your_module import Alignment, make_split_read_table
        >>> alignments = [
        ...     Alignment(read_name="read1", refname="chr1", start=100, end=150,
        ...               strand='+', clip1=10, match=40, clip2=5, pclip1=10),
        ...     Alignment(read_name="read2", refname="chr2", start=200, end=250,
        ...               strand='-', clip1=15, match=35, clip2=10, pclip1=15)
        ... ]
        >>> df = make_split_read_table(alignments)
        >>> print(df)
           qname chrom  start  end strand  clip1  match  clip2  pclip1
        0  read1  chr1    100  150      +     10     40      5      10
        1  read2  chr2    200  250      -     15     35     10      15
    """
    sa_cols = ['qname', 'chrom', 'start', 'end', 'strand', 'clip1', 'match', 'clip2', 'pclip1']
    data = []
    for alignment in alignments:
        qname = alignment.read_name
        chrom = alignment.refname
        start = alignment.start
        end = alignment.end
        strand = alignment.strand
        clip1 = alignment.clip1
        clip2 = alignment.clip2
        match = alignment.match
        pclip1 = alignment.pclip1
        field = [qname, chrom, start, end, strand, clip1, match, clip2, pclip1]
        data.append(field)
        # df.loc[df.shape[0]] = field
    df = pd.DataFrame(data, columns=sa_cols)
    df.drop_duplicates(inplace=True)
    df.sort_values(by=['qname', 'pclip1'], inplace=True)
    return df.reset_index(drop=True)

def is_breakpoints_not_sorted(chrom1, pos1, chrom2, pos2, chrom_order):
    if chrom_order.index(chrom1) > chrom_order.index(chrom2):
        return True
    if chrom_order.index(chrom1) == chrom_order.index(chrom2) and pos1 > pos2:
        return True
    return False

def reverse_complement(seq):
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'.'}
    revcmp = ''.join([comp[b] for b in seq[::-1]])
    return revcmp