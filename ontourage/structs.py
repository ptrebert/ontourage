from dataclasses import dataclass, field, astuple, asdict
import hashlib as hlib
import collections as col
import re


@dataclass
class Node:
    name: str
    length: int = 0
    orientation: int = 1
    sequence_md5: str = None
    # attributes related to graph topology
    # cc = connected component
    deg_in_plus: int = field(init=False, default=0)
    deg_in_minus: int = field(init=False, default=0)
    deg_out_plus: int = field(init=False, default=0)
    deg_out_minus: int = field(init=False, default=0)
    deg_total: int = field(init=False, default=0)
    cc_id: int = field(init=False, default=-1)
    cc_cardinality: int = field(init=False, default=0)
    centrality: float = field(init=False, default=-1)
    # attributes describing the node sequence composition
    num_A: int = field(init=False, default=0)
    pct_A: float = field(init=False, default=0)
    num_C: int = field(init=False, default=0)
    pct_C: float = field(init=False, default=0)
    num_G: int = field(init=False, default=0)
    pct_G: float = field(init=False, default=0)
    num_T: int = field(init=False, default=0)
    pct_T: float = field(init=False, default=0)
    num_N: int = field(init=False, default=0)
    pct_N: float = field(init=False, default=0)
    num_other: int = field(init=False, default=0)
    pct_other: float = field(init=False, default=0)
    num_bp_GvA: int = field(init=False, default=0)
    pct_bp_GvA: float = field(init=False, default=0)
    num_pattern_GpA: int = field(init=False, default=0)
    pct_pattern_GpA: float = field(init=False, default=0)
    num_bp_CvG: int = field(init=False, default=0)
    pct_bp_CvG: float = field(init=False, default=0)
    num_pattern_CpG: int = field(init=False, default=0)
    pct_pattern_CpG: float = field(init=False, default=0)

    def __post_init__(self):
        if self.orientation not in [-1, 1]:
            raise ValueError(f"Node orientation must be 1 or -1, and not {self.orientation}")
        if self.length < 0:
            raise ValueError(f"Node length cannot be negative: {self.length}")
        if self.sequence_md5 is not None and self.length == 0:
            raise ValueError(f"Sequence MD5 checksum is set, but node [sequence] length is 0: [MD5] {self.sequence_md5}")
        if self.sequence_md5 is None:
            self.sequence_md5 = 'n/a'  # avoid that None needs to be dumped as DataFrame/HDF

    def add_sequence_composition(self, sequence):
        self.sequence_md5 = hlib.md5(sequence.encode('ascii')).hexdigest()
        self.length = len(sequence)
        _upper_seq = sequence.upper()
        count_bases = col.Counter(_upper_seq)
        other = 0
        for base, count in count_bases.items():
            if base in ['A', 'C', 'G', 'T', 'N']:
                setattr(self, f'num_{base}', count)
                setattr(self, f'pct_{base}', round(count / self.length * 100, 2))
            else:
                other += count
        self.num_other = other
        self.pct_other = round(other / self.length * 100, 2)
        # record stats for common dinucleotides of interest (CG/CpG and GA/GpA [HiFi dropouts])
        self.num_bp_GvA = self.num_G + self.num_A
        self.pct_bp_GvA = round(self.num_bp_GvA / self.length * 100, 2)
        self._compute_pattern_stats(_upper_seq, 'GpA', re.compile('(GA)+'), 2)
        self.num_bp_CvG = self.num_C + self.num_G
        self.pct_bp_CvG = round(self.num_bp_CvG / self.length * 100, 2)
        self._compute_pattern_stats(_upper_seq, 'CpG', re.compile('(CG)+'), 2)
        return

    def _compute_pattern_stats(self, sequence, pattern_name, fixed_pattern, pattern_length=None):
        if pattern_length is None:
            pat_len = len(pattern_name)
        else:
            assert isinstance(pattern_length, int)
            pat_len = pattern_length
        num_bp = sum([m.end() - m.start() for m in fixed_pattern.finditer(sequence)])
        num_pattern = num_bp // pat_len
        setattr(self, f"num_pattern_{pattern_name}", num_pattern)
        # number of substrings of length k in string of length n: n - k + 1
        setattr(self, f"pct_pattern_{pattern_name}", round(num_pattern / (self.length - pat_len + 1) * 100, 2))
        return

    def set_node_degree(self, edge_counts):
        self.deg_out_plus = edge_counts[('out', 1)]
        self.deg_out_minus = edge_counts[('out', -1)]
        self.deg_in_plus = edge_counts[('in', 1)]
        self.deg_in_minus = edge_counts[('in', -1)]
        self.deg_total = self.deg_out_plus + self.deg_out_minus + self.deg_in_plus + self.deg_in_minus
        return

    def as_dict(self):
        return asdict(self)

    def as_tuple(self):
        return astuple(self)


@dataclass
class Edge:
    _id: str = field(init=False, default=None)
    node_a: str
    a_orientation: int
    node_b: str
    b_orientation: int
    length: int
    quality: int
    read_name: str = 'undefined'
    edge_source: str = 'assembler'
    edge_type: str = 'pair'
    error_type: str = 'undefined'

    def __post_init__(self):
        edge_id = ''.join(
            (str(x) for x in [self.node_a, self.a_orientation, self.node_b, self.b_orientation])
        )
        edge_hash = hlib.md5(edge_id.encode('ascii')).hexdigest()
        self._id = edge_hash
        if self.a_orientation not in [1, -1]:
            raise ValueError(f"Node A orientation must be 1 or -1, and not {self.a_orientation}")
        if self.b_orientation not in [1, -1]:
            raise ValueError(f"Node B orientation must be 1 or -1, and not {self.b_orientation}")
        if self.edge_type not in ['pair', 'loop']:
            raise ValueError(f"Type of edge must be 'pair' or 'loop', and not {self.edge_type}")

    def as_dict(self):
        return asdict(self)

    def as_tuple(self):
        return astuple(self)

    @staticmethod
    def make_edge(node_a, node_b, kwargs):
        quality = (node_a.support_quality + node_b.support_quality) // 2
        if node_a.support_source != node_b.support_source:
            # this should always be the case for the time being
            edge_source = f"{node_a.support_source}|{node_b.support_source}"
        else:
            edge_source = node_a.support_source
        full_kwargs = {
            'node_a': node_a.name,
            'a_orientation': node_a.orientation,
            'node_b': node_b.name,
            'b_orientation': node_b.orientation,
            'quality': quality,
            'edge_source': edge_source
        }
        full_kwargs.update(kwargs)
        edge = Edge(**kwargs)
        return edge


@dataclass(frozen=True)
class NodeSupport:
    name: str
    orientation: int
    length: int
    support_read: str
    support_quality: int
    support_start: int
    support_end: int
    support_bp: int
    support_pct: float
    support_source: str

    def as_dict(self):
        return asdict(self)

    def as_tuple(self):
        return astuple(self)

    def as_prefixed_record(self, prefix, read_align_start, read_align_end):

        record = dict((f'{prefix}_{key.replace("support", "align")}', val) for key, val in asdict(self).items())
        record[f'node_{prefix}'] = record['name']
        record[f'{prefix}_read_align_start'] = read_align_start
        record[f'{prefix}_read_align_end'] = read_align_end
        del record['name']
        return record


@dataclass(frozen=True)
class IssueRecord:
    node_a: str
    a_orientation: int
    a_length: int
    a_align_source: str  # alignment / GAF record ID
    a_align_read: str
    a_read_align_start: int
    a_read_align_end: int
    a_align_start: int
    a_align_end: int
    a_align_bp: int
    a_align_pct: float
    a_align_quality: int
    node_b: str
    b_orientation: int
    b_length: int
    b_align_source: str  # alignment / GAF record ID
    b_align_read: str
    b_read_align_start: int
    b_read_align_end: int
    b_align_start: int
    b_align_end: int
    b_align_bp: int
    b_align_pct: float
    b_align_quality: int
    description = str
    ovl_read_space_bp = int
    ovl_read_space_pct = float
    ovl_node_space_bp = int
    ovl_node_space_pct = float

    def as_dict(self):
        return asdict(self)

    def as_tuple(self):
        return astuple(self)

    @staticmethod
    def _get_overlap_keys():
        overlaps_keys = [
            'ovl_read_space_bp',
            'ovl_read_space_pct',
            'ovl_node_space_bp',
            'ovl_node_space_pct'
        ]
        return overlaps_keys

    @staticmethod
    def make_issue_record(record_a, record_b, description, overlaps):
        kwargs = {
            'description': description
        }
        if isinstance(overlaps, dict):
            kwargs.update(overlaps)
        else:
            ovl_keys = IssueRecord._get_overlap_keys()
            kwargs.update(dict((k, ovl) for k, ovl in zip(ovl_keys, overlaps)))
        kwargs.update(record_a)
        kwargs.update(record_b)
        issue = IssueRecord(**kwargs)
        return issue
