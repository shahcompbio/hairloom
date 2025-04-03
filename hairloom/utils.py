from typing import Set

import numpy as np
import pandas as pd

from hairloom.datatypes import Breakpoint, BreakpointChain

def enumerate_breakpoints(df):
    """Enumerates breakpoints from a DataFrame of genomic fragments.

    This function generates a ``BreakpointChain`` object containing ``Breakpoint`` objects
    based on the start and end positions of fragments in the input DataFrame. Breakpoints
    at the beginning and end of reads are treated differently to omit read boundaries.

    Args:
        df (pandas.DataFrame): A DataFrame containing fragment information, with the following columns:

            - `chrom`: The chromosome name (str).
            - `start`: The start position of the fragment (int).
            - `end`: The end position of the fragment (int).
            - `strand`: The strand information ('+' or '-').

    Returns:
        BreakpointChain: A chain of Breakpoint objects enumerated from the DataFrame.

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
    brks = []
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
    brks = BreakpointChain(brks)
    return brks

def get_secondaries(read):
    """Extracts secondary alignments from a sequencing read's 'SA' tag.

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
    """Creates a table summarizing split-read alignments.

    This function processes a list of alignment objects, extracting key attributes
    and organizing them into a structured pandas DataFrame. The resulting table
    is deduplicated and sorted for consistent organization.

    Args:
        alignments (list): A list of alignment objects. Each alignment object is expected
            to have the following attributes:

        - `read_name` (str): Query name of the read.
        - `refname` (str): Reference sequence name (chromosome or contig).
        - `start` (int): Start position of the alignment on the reference.
        - `end` (int): End position of the alignment on the reference.
        - `strand` (str): Strand information ('+' or '-').
        - `clip1` (int): Length of the first soft/hard clip before the matched region.
        - `clip2` (int): Length of the second soft/hard clip after the matched region.
        - `match` (int): Total number of matched bases in the alignment.
        - `pclip1` (int): Strand-corrected length of the first clip.

    Returns:
        pandas.DataFrame: A DataFrame containing the following columns:

        - `qname` (str): Query name of the read.
        - `chrom` (str): Reference sequence name (chromosome or contig).
        - `start` (int): Start position of the alignment on the reference.
        - `end` (int): End position of the alignment on the reference.
        - `strand` (str): Strand information ('+' or '-').
        - `clip1` (int): Length of the first soft/hard clip before the matched region.
        - `match` (int): Total number of matched bases in the alignment.
        - `clip2` (int): Length of the second soft/hard clip after the matched region.
        - `pclip1` (int): Strand-corrected length of the first clip.

    Notes:
        - Duplicate rows are removed based on the full set of columns.
        - The DataFrame is sorted by `qname` and `pclip1` for consistent organization.

    Example:
        >>> from hairloom.utils import make_split_read_table
        >>> class Alignment:
        ...     def __init__(self, read_name, refname, start, end, strand, clip1, match, clip2, pclip1):
        ...         self.read_name = read_name
        ...         self.refname = refname
        ...         self.start = start
        ...         self.end = end
        ...         self.strand = strand
        ...         self.clip1 = clip1
        ...         self.match = match
        ...         self.clip2 = clip2
        ...         self.pclip1 = pclip1
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


def make_seg_table(bundle, chroms=None):
    """Creates a table of genomic segments with support information.

    This function processes a bundle of ``BreakpointChain`` objects, extracts genomic segment
    coordinates, and combines them with support data into a structured pandas DataFrame.

    Args:
        bundle (list[BreakpointChain]): A list of ``BreakpointChain`` objects containing 
            segment information.
        chroms (list[str], optional): A list of chromosomes to include. If None, all chromosomes
            are considered. Defaults to None.

    Returns:
        pandas.DataFrame: A DataFrame containing the breakpoint table with the following
            columns:

        - `chrom` (str): Chromosome name.
        - `pos` (int): Position of the breakpoint.
        - `ori` (str): Orientation of the breakpoint ('+' or '-').
        - `support` (int): Support count for the breakpoint.

    Notes:
        - Segments are filtered based on the `chroms` parameter if provided.
        - Infinite values in the resulting DataFrame are replaced with NaN for consistency.

    Example:
        >>> from your_module import BreakpointChain, make_seg_table
        >>> bundle = [
        ...     Breakpoint('chr1', 100, '+'),
        ...     Breakpoint('chr2', 200, '+'),
        ...     Breakpoint('chr2', 300, '+'),
        ...     Breakpoint('chr1', 400, '+'),
        ... ]  # List of BreakpointChain objects
        >>> seg_df = make_seg_table(bundle, chroms=["chr1", "chr2"])
        >>> print(seg_df.head())
           chrom  pos1  pos2  support
        0  chr2   200   300       1
    """
    data = []
    
    for brks in bundle:
        for seg in brks.segs:
            assert seg.brk1.chrom == seg.brk2.chrom
            chrom = seg.brk1.chrom
            if chroms:
                if chrom not in chroms: 
                    continue
            pos1 = seg.brk1.pos
            pos2 = seg.brk2.pos
            ori1 = seg.brk1.ori
            ori2 = seg.brk2.ori
            if pos1 > pos2:
                pos1, ori1, pos2, ori2 = pos2, ori2, pos1, ori1
            # assert ori1 == '-' and ori2 == '+', (ori1, ori2)
            field = [chrom, pos1, pos2]
            data.append(field)
    seg_cols = ['chrom', 'pos1', 'pos2']
    seg_df = pd.DataFrame(data, columns=seg_cols)
    seg_df = seg_df.value_counts().reset_index().rename(columns={'count':'support', 0:'support'})
    seg_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    return seg_df


def make_brk_table(bundle, chroms=None):
    """Generates a table of breakpoints with support information.

    This function processes a bundle of ``BreakpointChain`` objects, extracts
    breakpoint coordinates, and organizes them into a structured pandas
    DataFrame. Optionally, the resulting table can be filtered by specific
    chromosomes.

    Args:
        bundle (list[BreakpointChain]): A list of ``BreakpointChain`` objects containing
            breakpoint information.
        chroms (list[str], optional): A list of chromosome names to include. If None, all
            chromosomes are considered. Defaults to None.

    Returns:
        pandas.DataFrame: A DataFrame containing the breakpoint table with the following
            columns:

        - `chrom` (str): Chromosome name.
        - `pos` (int): Position of the breakpoint.
        - `ori` (str): Orientation of the breakpoint ('+' or '-').
        - `support` (int): Support count for the breakpoint.

    Notes:
        - Breakpoints are filtered based on the `chroms` parameter if provided.
        - Infinite values in the resulting DataFrame are replaced with NaN to
          maintain consistency.

    Example:
        >>> from hairloom.datatypes import Breakpoint, BreakpointChain
        >>> from hairloom.utils import make_brk_table
        >>> bundle = [
        ...     BreakpointChain([
        ...         Breakpoint("chr1", 100, "+"),
        ...         Breakpoint("chr1", 200, "+"),
        ...         Breakpoint("chr2", 300, "+"),
        ...     ]),
        ... ]
        >>> chroms = ["chr1", "chr2"]
        >>> brk_df = make_brk_table(bundle, chroms)
        >>> print(brk_df)
           chrom  pos ori  support
        0  chr1  100   +       1
        1  chr1  200   +       1
        2  chr2  300   +       1
    """
    data = []
    for brks in bundle:
        for brk in brks:
            if chroms:
                if brk.chrom not in chroms:
                    continue
            field = (brk.chrom, brk.pos, brk.ori)
            data.append(field)

    brk_cols = ['chrom', 'pos', 'ori']
    brk_df = pd.DataFrame(data, columns=brk_cols)
    brk_df = brk_df.value_counts().reset_index().rename(columns={'count':'support', 0:'support'})
    brk_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    
    return brk_df


def make_tra_table(bundle):
    """Creates a table of translocations with support information.

    This function processes a bundle of ``BreakpointChain`` objects, extracts translocation
    coordinates, and combines them with support data into a structured pandas DataFrame.

    Args:
        bundle (list[BreakpointChain]): A list of ``BreakpointChain`` objects containing 
            translocation information.

    Returns:
        pandas.DataFrame: A DataFrame containing the following columns:

        - `chrom1` (str): Chromosome name of the first breakpoint.
        - `pos1` (int): Position of the first breakpoint.
        - `ori1` (str): Orientation of the first breakpoint ('+' or '-').
        - `chrom2` (str): Chromosome name of the second breakpoint.
        - `pos2` (int): Position of the second breakpoint.
        - `ori2` (str): Orientation of the second breakpoint ('+' or '-').
        - `support` (int): Support count for the translocation.

    Notes:
        - Translocations are identified by pairs of breakpoints (`BreakpointPair` objects) in
          the `tras` attribute of each `BreakpointChain`.
        - Duplicate translocations are removed by maintaining a `set` of seen coordinate pairs.
        - Infinite values in the resulting DataFrame are replaced with NaN for consistency.

    Example:
        >>> from hairloom.datatypes import Breakpoint, BreakpointChain
        >>> from hairloom.utils import make_tra_table
        >>> chain1 = BreakpointChain([
        ...     Breakpoint('chr1', 100, '+'),
        ...     Breakpoint('chr2', 200, '+')
        ... ])
        >>> chain2 = BreakpointChain([
        ...     Breakpoint('chr2', 300, '+'),
        ...     Breakpoint('chr1', 400, '+')
        ... ])
        >>> bundle = [chain1, chain2]
        >>> tra_df = make_tra_table(bundle)
        >>> print(tra_df)
           chrom1  pos1 ori1 chrom2  pos2 ori2  support
        0   chr1   100    +   chr2   200    +       1
        1   chr2   300    +   chr1   400    +       1
    """
    data = []
    for brks in bundle:
        for tra in brks.tras:
            coord1 = (tra.brk1.chrom, tra.brk1.pos, tra.brk1.ori)
            coord2 = (tra.brk2.chrom, tra.brk2.pos, tra.brk2.ori)
            field = [*coord1, *coord2]
            data.append(field)
    
    tra_cols = ['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2']
    tra_df = pd.DataFrame(data, columns=tra_cols)
    tra_df = tra_df.value_counts().reset_index().rename(columns={'count':'support', 0:'support'})
    tra_df.replace([-np.inf, np.inf], np.nan, inplace=True) 
    return tra_df


def melt_tra_table(df: pd.DataFrame) -> Set[Breakpoint]:
    """Converts a DataFrame of structural variants (SVs) into a set of unique breakpoints.

    This function processes a DataFrame with columns representing two breakpoints per row 
    (`chrom1`, `pos1`, `ori1`, `chrom2`, `pos2`, `ori2`), extracts these breakpoints, and 
    returns them as a set of `Breakpoint` objects.

    Args:
        df (pandas.DataFrame): A DataFrame containing structural variant information 
            with the following columns:
            - `chrom1` (str): Chromosome of the first breakpoint.
            - `pos1` (int): Position of the first breakpoint.
            - `ori1` (str): Orientation of the first breakpoint ('+' or '-').
            - `chrom2` (str): Chromosome of the second breakpoint.
            - `pos2` (int): Position of the second breakpoint.
            - `ori2` (str): Orientation of the second breakpoint ('+' or '-').

    Returns:
        set[Breakpoint]: A set of unique `Breakpoint` objects representing all breakpoints 
        found in the DataFrame.

    Notes:
        - Each row in the DataFrame corresponds to a pair of breakpoints, and this function
          extracts both breakpoints from each row.
        - Duplicate breakpoints across rows are automatically removed by the set.

    Example:
        >>> import pandas as pd
        >>> from hairloom.datatypes import Breakpoint
        >>> from hairloom.utils import melt_tra_table
        >>> data = {
        ...     'chrom1': ['chr1', 'chr2'],
        ...     'pos1': [100, 200],
        ...     'ori1': ['+', '-'],
        ...     'chrom2': ['chr3', 'chr4'],
        ...     'pos2': [300, 400],
        ...     'ori2': ['-', '+']
        ... }
        >>> df = pd.DataFrame(data)
        >>> result = melt_tra_table(df)
        >>> print(result)
        {Breakpoint('chr1', 100, '+'), Breakpoint('chr3', 300, '-'),
         Breakpoint('chr2', 200, '-'), Breakpoint('chr4', 400, '+')}
    """
    breakpoints = set(
        Breakpoint(row['chrom1'], row['pos1'], row['ori1']) 
        for row in df[['chrom1', 'pos1', 'ori1']].to_dict('records')
    ).union(
        Breakpoint(row['chrom2'], row['pos2'], row['ori2']) 
        for row in df[['chrom2', 'pos2', 'ori2']].to_dict('records')
    )
    return breakpoints


def melt_brk_table(df: pd.DataFrame) -> Set[Breakpoint]:
    """
    Converts a DataFrame of breakpoints into a set of unique `Breakpoint` objects.

    This function processes a DataFrame containing breakpoint data (chromosome, position, 
    orientation, and support) and returns a set of unique `Breakpoint` objects representing 
    all breakpoints in the table.

    Args:
        df (pandas.DataFrame): A DataFrame containing breakpoint information with the 
            following columns:
            - `chrom` (str): Chromosome name of the breakpoint.
            - `pos` (int): Position of the breakpoint.
            - `ori` (str): Orientation of the breakpoint ('+' or '-').
            - `support` (int): Support count for the breakpoint.

    Returns:
        set[Breakpoint]: A set of unique `Breakpoint` objects representing all breakpoints
        in the DataFrame.

    Notes:
        - The `support` column is ignored in the resulting `Breakpoint` set, as it is 
          not part of the `Breakpoint` attributes.
        - Duplicate breakpoints in the DataFrame are automatically removed by the set.

    Example:
        >>> import pandas as pd
        >>> from hairloom.datatypes import Breakpoint
        >>> from hairloom.utils import melt_brk_table
        >>> data = {
        ...     'chrom': ['chr1', 'chr1', 'chr2'],
        ...     'pos': [10, 100, 20],
        ...     'ori': ['+', '-', '+'],
        ...     'support': [10, 5, 2],
        ... }
        >>> df = pd.DataFrame(data)
        >>> result = melt_brk_table(df)
        >>> print(result)
        {Breakpoint('chr1', 10, '+'), Breakpoint('chr1', 100, '-'), Breakpoint('chr2', 20, '+')}
    """
    # Use pandas vectorized operations to construct Breakpoint objects
    return {
        Breakpoint(chrom, pos, ori)
        for chrom, pos, ori in zip(df['chrom'], df['pos'], df['ori'])
    }