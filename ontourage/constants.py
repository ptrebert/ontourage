
LOG_MESSAGE_FORMAT = "%(asctime)s | %(funcName)s | %(levelname)s | %(message)s"

ORIENTATION_MAP = {
    '>': 1,
    '<': -1,
    '+': 1,
    '-': -1,
}

REVERSE_ORIENTATION_MAP = {
    1: '+',
    -1: '-'
}

REVCOMP_MAPPING = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}
REVCOMP_MAPPING.update({k.lower(): v.lower() for k, v in REVCOMP_MAPPING.items()})

REVCOMP_TABLE = str.maketrans(REVCOMP_MAPPING)

SEGMENT_ORIENTATION_SYMBOLS = '(<|>)'

GAF_DEFAULT_COLUMNS = [
    'read_name',
    'read_length',
    'read_align_start',
    'read_align_end',
    'read_align_orientation',
    'path',
    'path_length',
    'path_align_start',
    'path_align_end',
    'num_matches',
    'align_block_length',
    'mapq',
]

GAF_DEFAULT_COLUMN_TYPES = [
    str,
    int,
    int,
    int,
    None,
    str,
    int,
    int,
    int,
    int,
    int,
    int
]

GAF_DEFAULT_COLUMN_TYPE_MAP = dict((c, t) for c, t in zip(GAF_DEFAULT_COLUMNS, GAF_DEFAULT_COLUMN_TYPES) if t is not None)
