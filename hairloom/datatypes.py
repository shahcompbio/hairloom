import re
import pyfaidx


class SplitAlignment:
    """Represents a split alignment from a sequencing read.

    Parses the CIGAR string and extracts alignment information such as clip lengths,
    matched bases, and strand-corrected values.

    Attributes:
        read_name (str): The name of the sequencing read.
        refname (str): The reference sequence name (e.g., chromosome or contig).
        cigarstring (str): The CIGAR string of the alignment.
        start (int): The start position of the alignment on the reference.
        strand (str): The strand information ('+' or '-').

        cigar_tuples (list[tuple]): Parsed CIGAR string as a list of (operation, length) tuples.
        primary (NoneType): Placeholder for primary alignment information, initialized to `None`.
        match (int): Total number of matched bases in the alignment.
        aln_cols (list[str]): Column headers for alignment fields.
        clip1 (int): Length of the first clip (soft/hard) before the matched region.
        clip2 (int): Length of the second clip (soft/hard) after the matched region.
        end (int): The end position of the alignment on the reference.
        pclip1 (int): Strand-corrected length of the first clip.
    """
    def __init__(self, cigarstring, read_name, refname, read_pos, strand):
        """Initializes the SplitAlignment instance.

        Args:
            cigarstring (str): The CIGAR string of the alignment.
            read_name (str): The name of the sequencing read.
            refname (str): The reference sequence name (chromosome or contig).
            read_pos (int): The start position of the alignment on the reference.
            strand (str): The strand information ('+' or '-').

        Raises:
            ValueError: If the strand is not '+' or '-'.
        """
        self.read_name = read_name
        self.refname = refname
        self.cigarstring = cigarstring
        self.cigar_tuples = SplitAlignment.get_cigar_tuples(self.cigarstring)
        self.primary = None
        self.start = read_pos
        self.strand = strand
        self.match = 0

        self.aln_cols = ['qname', 'chrom', 'start', 'end', 'strand', 'clip1', 'match', 'clip2', 'pclip1']
        self.extract_cigar_field()
        
    def extract_cigar_field(self):
        """Parses the CIGAR string to calculate clip lengths, matched bases, and alignment end position."""
        split_ctags = {4, 5}
        query_consumers = {0, 1, 4, 7, 8}
        reference_consumers = {0, 2, 3, 7, 8}
        self.clip1, self.match, self.clip2 = 0, 0, 0
        flag_after_match = False
        chrom, pos, strand, cigar_tuples = self.refname, self.start, self.strand, self.cigar_tuples
        tpos = int(pos) # init reference pos
        for ctag, clen in cigar_tuples:
            if ctag in split_ctags:
                if not flag_after_match:
                    self.clip1 += clen
                else:
                    self.clip2 += clen
            elif ctag in query_consumers:
                flag_after_match = True
                self.match += clen
            if ctag in reference_consumers:
                flag_after_match = True
                tpos += clen
                continue
        self.end = tpos
        pclip1 = self.clip1 if strand == '+' else self.clip2
        self.pclip1 = pclip1 # positively corrected first clip length

    @staticmethod
    def get_cigar_tuples(cigarstring):
        """Parses a CIGAR string and converts it into a list of operation-length tuples.

        Args:
            cigarstring (str): The CIGAR string to parse, following the standard format
                used in sequence alignments (e.g., "10M5I20M").

        Returns:
            list[tuple[int, int]]: Parsed CIGAR operations and lengths.

        Raises:
            ValueError: If the CIGAR string contains invalid operations.

        Example:
            >>> cigarstring = "10M5I20M"
            >>> SplitAlignment.get_cigar_tuples(cigarstring)
            [(0, 10), (1, 5), (0, 20)]

        Notes:
            Supported CIGAR operations:
                - 'M' (0): Alignment match (can be a sequence match or mismatch).
                - 'I' (1): Insertion to the reference.
                - 'D' (2): Deletion from the reference.
                - 'N' (3): Skipped region from the reference.
                - 'S' (4): Soft clipping (clipped sequences present in the read).
                - 'H' (5): Hard clipping (clipped sequences not present in the read).
                - 'P' (6): Padding (silent deletion from padded reference).
                - '=' (7): Sequence match.
                - 'X' (8): Sequence mismatch.
        """
        cigar_types = 'MIDNSHP=X'
        converter = {
            cigar_types[i]:i for i in range(len(cigar_types))
        }
        cigar_tuples = []
        for match in re.finditer(r'(\d+)([A-Z=])', cigarstring):
            clen, ctype = match.groups()
            cigar_tuple = (converter[ctype], int(clen))
            cigar_tuples.append(cigar_tuple)
        return cigar_tuples


class Breakpoint:
    """Represents a genomic breakpoint with associated properties and methods.

    Attributes:
        chrom (str): The chromosome name where the breakpoint is located.
        pos (int): The 1-based position of the breakpoint on the chromosome.
        ori (str): The orientation of the breakpoint ('+' or '-').

        upstream (str or None): Sequence upstream of the breakpoint, initialized to None.
        downstream (str or None): Sequence downstream of the breakpoint, initialized to None.
        seq_rearranged (str or None): Rearranged sequence at the breakpoint, initialized to None.
        seq_removed (str or None): Removed sequence at the breakpoint, initialized to None.

        chroms (list[str]): List of valid chromosome names, including both standard
            ('1', '2', ..., 'X', 'Y', 'M') and prefixed ('chr1', 'chr2', ..., 'chrX', 'chrY').
    """
    chroms = [str(c) for c in range(1, 22+1)] + ['X', 'Y', 'M']
    chroms += ['chr'+c for c in chroms]
    def __init__(self, chrom, pos, orientation):
        """Initializes a Breakpoint instance.

        Args:
            chrom (str): The chromosome name where the breakpoint is located.
            pos (int): The 1-based position of the breakpoint on the chromosome.
            ori (str): The orientation of the breakpoint ('+' or '-').

        Raises:
            ValueError: If `ori` is not '+' or '-'.
        """
        self.chrom = chrom
        self.pos = pos
        self.ori = orientation
        self.upstream = None
        self.downstream = None
        self.seq_rearranged = None
        self.seq_removed = None

    def get_breakpoint_seqs(self, margin:int, genome:pyfaidx.Fasta):
        """Retrieves upstream and downstream sequences around the breakpoint.

        Computes rearranged and removed sequences based on the given margin
        and genome dictionary.

        Args:
            margin (int): Number of bases to include upstream and downstream of the breakpoint.
            genome (dict): A dictionary mapping chromosome names to their respective sequences.

        Raises:
            ValueError: If the chromosome is not found in the genome.
        """
        self.upstream, self.downstream = get_breakpoint_seqs(self.chrom, self.pos, margin, genome)
        if self.ori == '+':
            self.seq_rearranged = self.upstream
            self.seq_removed = self.downstream
        elif self.ori == '-':
            self.seq_rearranged = self.downstream
            self.seq_removed = self.upstream
        else:
            raise ValueError(f'self.ori = {self.ori}')

    def __repr__(self):
        return f'{self.chrom}:{self.pos}:{self.ori}'

    def __lt__(self, other):
        self_chrom_ix, other_chrom_ix = None, None
        if self.chrom in self.chroms:
            self_chrom_ix = self.chroms.index(self.chrom)
        if other.chrom in self.chroms:
            other_chrom_ix = self.chroms.index(other.chrom)
        if self_chrom_ix and not other_chrom_ix:
            return False
        if not self_chrom_ix and other_chrom_ix:
            return True
        if not self_chrom_ix and not other_chrom_ix:
            return self.chrom < other.chrom
        if self_chrom_ix < other_chrom_ix: # only if both not None
            return True
        elif self_chrom_ix == other_chrom_ix and self.pos < other.pos:
            return True
        return False
    
    def __eq__(self, other):
        if not isinstance(other, Breakpoint):
            return NotImplemented
        return (
            self.chrom == other.chrom and
            self.pos == other.pos and
            self.ori == other.ori
        )

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ori))


class BreakpointPair:
    """Represents a pair of genomic breakpoints.

    Attributes:
        brk1 (Breakpoint): The first breakpoint in the pair.
        brk2 (Breakpoint): The second breakpoint in the pair.

        aln_segment (bool): Indicates whether this pair is part of an alignment segment.
            Defaults to `False`.

    Methods:
        __repr__(): Returns a string representation of the breakpoint pair.
    """
    def __init__(self, brk1, brk2):
        self.brk1 = brk1
        self.brk2 = brk2
        self.aln_segment = False

    def __repr__(self):
        return f'{self.brk1.chrom}:{self.brk1.pos}:{self.brk1.ori}-{self.brk2.chrom}:{self.brk2.pos}:{self.brk2.ori}'

    def __eq__(self, other):
        if not isinstance(other, BreakpointPair):
            return NotImplemented
        return (
            (self.brk1 == other.brk1 and self.brk2 == other.brk2) or
            (self.brk1 == other.brk2 and self.brk2 == other.brk1)  # Allow unordered equality
        )

    def __hash__(self):
        return hash(frozenset({self.brk1, self.brk2}))
    

class BreakpointChain(list):
    """Represents a chain of genomic breakpoints.

    This will often represent the chain of breakpoints coming from a single read.
    This class extends the Python list to store breakpoints and provides methods
    to enumerate transitions and segments.

    Attributes:
        tras (list[BreakpointPair]): List of transitions (pairs of breakpoints).
        segs (list[BreakpointPair]): List of segments (pairs of breakpoints).

    Args:
        brks_iterable (iterable): An iterable containing breakpoint objects.

    Example:
        >>> brk1 = Breakpoint("chr1", 100, "+")
        >>> brk2 = Breakpoint("chr1", 200, "-")
        >>> chain = BreakpointChain([brk1, brk2])
        >>> chain.tras
        [BreakpointPair(brk1, brk2)]
    """
    def __init__(self, brks_iterable):
        if len(brks_iterable) % 2 != 0:
            raise ValueError("BreakpointChain length must be a multiple of two.")
        super().__init__(brks_iterable)
        self._get_transitions()
        self._get_segments()
        self.qname = None # read qname slot

    # enumeration of tras
    def _get_transitions(self, sort_transition=True):
        """Enumerates transitions in the chain.

        Args:
            sort_transition (bool): If True, ensures transitions are sorted
                such that the smaller breakpoint comes first. Defaults to True.
        """
        self.tras = []
        if len(self) >= 2:
            ix_range = range(0, len(self), 2)
            for i in ix_range:
                brk1 = self[i]
                brk2 = self[i+1]
                if sort_transition:
                    if brk1 > brk2:
                        brk1, brk2 = brk2, brk1
                tra = BreakpointPair(brk1, brk2)
                self.tras.append(tra)

    # enumeration of segs
    def _get_segments(self):
        """
        Enumerates segments (pairs of breakpoints) in the chain and appends
        them to the `segs` attribute.

        Notes:
            - Segments are created from consecutive breakpoints in the chain,
              skipping the first and last breakpoints.
            - If the chain has fewer than 3 breakpoints, no segments are added.

        Example:
            >>> chain = BreakpointChain([brk1, brk2, brk3, brk4, brk5])
            >>> chain.segs
            [BreakpointPair(brk2, brk3), BreakpointPair(brk4, brk5)]
        """
        self.segs = []
        ix_range = range(1, len(self)-1, 2)
        for i in ix_range:
            brk1 = self[i]
            brk2 = self[i+1]
            if brk1 > brk2:
                brk1, brk2 = brk2, brk1
            # assert brk1.ori == '-' and brk2.ori == '+' # --- [ READ SEQUENCE ] --- form
            seg = BreakpointPair(brk1, brk2)
            self.segs.append(seg)


class Transitions:
    """Calculates transitions (tra) between genomic fragments.

    Attributes:
        list (list[tuple]): A list of transitions, where each transition
            is a tuple of the form ((chrom1, pos1, ori1), (chrom2, pos2, ori2)).

    Args:
        df (pandas.DataFrame): A DataFrame with fragment information,
            including 'qname', 'chrom', 'start', 'end', and 'strand'.

    Example:
        >>> df = pd.DataFrame({
        ...     'qname': ['read1', 'read1'],
        ...     'chrom': ['chr1', 'chr2'],
        ...     'start': [100, 200],
        ...     'end': [150, 250],
        ...     'strand': ['+', '-']
        ... })
        >>> transitions = Transitions(df)
        >>> transitions.list
        [(('chr1', 150, '+'), ('chr2', 200, '-'))]
    """
    def __init__(self, df):
        self.df = df.copy()
        self.list = []
        self.get_list()

    def get_list(self):
        """
        Computes transitions (tra) from the DataFrame.

        Notes:
            - Transitions are calculated between consecutive fragments in the DataFrame 
              grouped by 'qname'.
            - For each fragment pair, the orientation and positions are determined based 
              on the strand, creating a transition tuple of the form 
              ((chrom1, pos1, ori1), (chrom2, pos2, ori2)).

        Modifies:
            list (list[tuple]): Appends computed transitions to the `list` attribute.

        Example:
            >>> df = pd.DataFrame({
            ...     'qname': ['read1', 'read1', 'read2'],
            ...     'chrom': ['chr1', 'chr2', 'chr3'],
            ...     'start': [100, 200, 300],
            ...     'end': [150, 250, 350],
            ...     'strand': ['+', '-', '+']
            ... })
            >>> transitions = Transitions(df)
            >>> transitions.list # Note: 'read2' won't be included
            [(('chr1', 150, '+'), ('chr2', 250, '-'))]
        """
        for qname, qdf in self.df.groupby('qname'):
            qdf.reset_index(drop=True, inplace=True)
            n_fragments = qdf.shape[0]
            if n_fragments <= 1: continue # don't count TRA when none
            prevdf = qdf.shift(1)
            for rix, row in qdf.iterrows():
                if rix == 0: continue
                prow = prevdf.loc[rix]
                chrom1, start1, end1, strand1 = prow['chrom'], int(prow['start']), int(prow['end']), prow['strand']
                chrom2, start2, end2, strand2 = row['chrom'], int(row['start']), int(row['end']), row['strand']
                if strand1 == '+':
                    ori1 = '+'
                    pos1 = end1
                else:
                    ori1 = '-'
                    pos1 = start1
                if strand2 == '+':
                    ori2 = '-'
                    pos2 = start2
                else:
                    ori2 = '+'
                    pos2 = end2
                tra = ((chrom1, pos1, ori1), (chrom2, pos2, ori2))
                self.list.append(tra)
                

class Segments:
    """Calculates and stores genomic segments (middle fragments).

    Attributes:
        list (list[tuple]): A list of segments, each represented as a tuple
            (chrom, start, end).

    Args:
        df (pandas.DataFrame): A DataFrame with fragment information,
            including 'qname', 'chrom', 'start', and 'end'.

    Example:
        >>> df = pd.DataFrame({
        ...     'qname': ['read1', 'read1', 'read1'],
        ...     'chrom': ['chr1', 'chr1', 'chr1'],
        ...     'start': [100, 200, 300],
        ...     'end': [150, 250, 350],
        ... })
        >>> segments = Segments(df)
        >>> segments.list
        [('chr1', 200, 250)]
    """
    def __init__(self, df):
        self.df = df.copy()
        self.list = []
        self.get_list()

    def get_list(self):
        """
        Computes segments (middle fragments) from the DataFrame.

        Notes:
            - Segments are defined as genomic fragments that are neither the first 
              nor the last fragment in a group (grouped by 'qname').
            - If a group contains fewer than three fragments, no segments are added.

        Modifies:
            list (list[tuple]): Appends computed segments to the `list` attribute.

        Example:
            >>> df = pd.DataFrame({
            ...     'qname': ['read1', 'read1', 'read1', 'read2'],
            ...     'chrom': ['chr1', 'chr1', 'chr1', 'chr2'],
            ...     'start': [100, 200, 300, 400],
            ...     'end': [150, 250, 350, 450],
            ... })
            >>> segments = Segments(df)
            >>> segments.list
            [('chr1', 200, 250)]
        """
        for qname, qdf in self.df.groupby('qname'):
            n_fragments = qdf.shape[0]
            if n_fragments <= 2: continue # don't count segments when none
            for rix, row in qdf.iterrows():
                if rix == 0 or rix == n_fragments-1:
                    continue
                chrom, start, end = row['chrom'], int(row['start']), int(row['end'])
                segment = (chrom, start, end)
                self.list.append(segment)


def get_breakpoint_seqs(chrom, pos, margin, genome):
    """Extracts upstream and downstream sequences around a breakpoint.

    Args:
        chrom (str): Chromosome name.
        pos (int): 1-based position of the breakpoint.
        margin (int): Number of bases upstream and downstream to extract.
        genome (dict): Dictionary mapping chromosome names to sequences.

    Returns:
        tuple[str, str]: A tuple containing:
            - `upstream`: Sequence upstream of the breakpoint.
            - `downstream`: Sequence downstream of the breakpoint.

    Example:
        >>> genome = {'chr1': "A" * 1000, 'chr2': "T" * 1000}
        >>> get_breakpoint_seqs('chr1', 5, 3, genome)
        ('AA', 'AAAA')
    """
    linear_chroms = Breakpoint.chroms
    assert chrom in genome, chrom
    zpos = pos-1 # zero-based
    chrom_len = len(genome[chrom])
    start = zpos - margin
    end = zpos + margin
    if chrom in linear_chroms:
        start = max(0, start)
        end = min(chrom_len, end)
    
    if start < 0: # for circular chromosomes, take `3'` as upstream
        upstream = genome[chrom][start:].seq.upper() + genome[chrom][:zpos].seq.upper() # upstream of 0
    else:
        upstream = genome[chrom][start : zpos].seq.upper() # upstream
    if end >= chrom_len: # for circular chromosomes, take `5'` as downstream
        residual = end - chrom_len
        downstream = genome[chrom][zpos:].seq.upper() + genome[chrom][:residual].seq.upper() # downstream over end
    else:
        downstream = genome[chrom][zpos : end].seq.upper() # downstream
    return upstream, downstream
