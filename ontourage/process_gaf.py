import collections as col
import argparse as argp
import pathlib as pl
import functools as fnt
import logging as logging
import itertools as itt
import multiprocessing as mp
import re
import pickle as pck

import pandas as pd
import xopen

from ontourage import GAF_DEFAULT_COLUMN_TYPE_MAP, ORIENTATION_MAP, GAF_DEFAULT_COLUMNS, GRAPH_SEGMENT_REGEXP

logger = logging.getLogger()


def process_gaf_cli_parser(module_parsers):

    name = 'align'
    desc = 'Process a graph alignment (GAF) input file and cache its data.'

    if module_parsers is None:
        parser = argp.ArgumentParser(prog=name, description=desc, add_help=True)
    else:
        parser = module_parsers.add_parser(
            name,
            description=desc,
            help=desc,
            add_help=True
        )

    io_group = parser.add_argument_group('Input/Output')

    io_group.add_argument(
        "--align",
        "--gaf",
        "-a",
        required=True,
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="align",
        metavar="GAF",
        help="Path to input GAF file. Can be compressed (gzip, bzip, xz)."
    )

    io_group.add_argument(
        "--cache-folder",
        "-cf",
        required=True,
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="cache_folder",
        metavar="CACHE-FOLDER",
        help="Path to an existing folder for caching. "
             "Note that the folder has to the same as used for the graph (GFA) caching."
    )
    io_group.add_argument(
        "--cache-prefix",
        "-cp",
        default="ONTourage",
        type=str,
        dest="cache_prefix",
        metavar="CACHE-PREFIX",
        help="Specify a prefix for the individual cache files. Default: ONTourage "
             "Note that the prefix has to match the prefix specified for the graph (GFA) caching."
    )
    io_group.add_argument(
        "--ignore-cache",
        "-i",
        action="store_true",
        default=False,
        dest="ignore_cache",
        help="Ignore existing cached data and rerun everything. "
             "Note that the graph (GFA) cached data are still required."
    )
    parser.set_defaults(execute=run_process_gaf)
    if module_parsers is None:
        return parser
    else:
        return module_parsers


def compute_supported_node_sequence(start_offset, node_length, edge_length, aln_start_in_path, aln_end_in_path):
    """
    Sequentially compute sequence positions of nodes in paths that are supported by a read alignment

    Example:
    - path w/ 3 nodes: 0 --- 100 ; 75 --- 175 ; 150 --- 250 = path length in bp is 250
    --- (= 3 calls to this function)
    - 1 alignment: 30 --- 230 = alignment length in bp is 200 (in path space)
    ---> first 30 bp in node 1 are not supported
    ---> last 20 bp in node 3 are not supported

    case node 1:
    - offset: 0 (first node in path)
    - node_start_in_path: 0 (offset)
    - node_end_in_path: 100 (offset + node_length)
    - supported sequence: positions 30 - 100, compute as follows:
    --- max(node_start_in_path, aln_start_in_path) - node_start_in_path
    --- max(0, 30) - 0 = 30
    --- (subtract node_start_in_path to shift to zero / local segment coordinates)
    --- min(node_end_in_path, aln_end_in_path) - node_start_in_path
    --- min(100, 230) - 0 = 100
    - result: positions 30 - 100 = 70 bp supported
    - adjust offset:
    --- node_end_in_path - edge_length = new offset for next node
    --- 100 - 25 = 75

    case node 2:
    - offset: 75
    - node_start_in_path: 75 (offset)
    - node_end_in_path: 175 (offset + node_length)
    - supported_sequence:
    --- max(75, 30) - 75 = 0
    --- min(175, 230) - 75 = 100
    - result: positions 0 - 100 = 100 bp supported
    adjust offset:
    - 175 - 25 = 150 (new offset for next node)

    case node 3:
    - offset: 150
    - node_start_in_path: 150
    - node_end_in_path: 150 + 100 = 250
    - supported sequence:
    --- max(150, 30) - 150 = 0
    --- min(250, 230) - 150 = 80
    - result: positions 0 - 80 = 80 bp supported

    """

    node_start_in_path = start_offset
    node_end_in_path = node_start_in_path + node_length

    node_support_start = max(node_start_in_path, aln_start_in_path) - node_start_in_path
    node_support_end = min(node_end_in_path, aln_end_in_path) - node_start_in_path
    assert node_support_end > node_support_start
    node_support_length = node_support_end - node_support_start
    node_support_pct = round(node_support_length / node_length * 100, 2)

    new_offset = node_end_in_path - edge_length
    align_stats = {
        'support_start': node_support_start,
        'support_end': node_support_end,
        'support_length': node_support_length,
        'support_pct': node_support_pct,
    }
    return align_stats, new_offset


def enumerate_overlapping_query_blocks(alignments):
    """
    Enumerate all alignment blocks relative to the query (= the read)
    to simplify identifying overlapping alignments.
    """
    if alignments.shape[0] == 1:
        alignments['read_aln_block'] = 0
    else:
        enum_blocks = alignments[['read_align_start', 'read_align_end']].sort_values('read_align_start', inplace=False)
        enum_blocks['read_aln_block'] = (enum_blocks['read_align_start'] > enum_blocks['read_align_end'].shift().cummax()).cumsum()
        alignments.loc[enum_blocks.index, 'read_aln_block'] = enum_blocks['read_aln_block']
    return alignments


def process_path_alignments(read_name, node_lengths, edge_lengths, segment_re, alignments):
    """
    Process all path alignments |P| > 1 for a single read
    - track path members
    - transform single alignment to series of node-pair alignments (= supported edges)
    - compute supported sequence per node in path
    """
    supported_edges = []
    supported_sequences = []
    path_nodes = set()
    nodes_per_block = col.defaultdict(set)

    match_segment = re.compile(segment_re)
    for _, row in alignments.iterrows():
        path_members = [m.group(0) for m in re.finditer(match_segment, row['path'])]
        start_offset = 0
        # in Python 3.10+: itt.pairwise(path_members):
        for (node_a, node_b) in zip(path_members[:-1], path_members[1:]):
            a_orient = ORIENTATION_MAP[node_a[0]]
            node_a = node_a[1:]
            path_nodes.add(node_a)
            nodes_per_block[row['read_aln_block']].add((node_a, row['mapq']))
            a_length = node_lengths[node_a]

            b_orient = ORIENTATION_MAP[node_b[0]]
            node_b = node_b[1:]
            path_nodes.add(node_b)
            nodes_per_block[row['read_aln_block']].add((node_b, row['mapq']))
            b_length = node_lengths[node_b]

            # alignment to path implies that edge exists in graph
            edge_length = edge_lengths[(node_a, a_orient, node_b, b_orient)]
            attributes = {
                'a_length': a_length,
                'a_orient': a_orient,
                'b_length': b_length,
                'b_orient': b_orient,
                'overlap': edge_length,
                'source': 'alignment',
                'quality': row['mapq'],
                'read_name': read_name,
            }
            supported_edges.append((node_a, node_b, attributes))
            if start_offset == 0:
                a_support, start_offset = compute_supported_node_sequence(
                    start_offset,
                    a_length,
                    edge_length,
                    row['path_align_start'],
                    row['path_align_end']
                )
                a_support['name'] = node_a
                a_support['orientation'] = a_orient
                a_support['read_name'] = read_name
                a_support['quality'] = row['mapq']
                a_support['length'] = a_length
                a_support['align_type'] = 'path_member'
                supported_sequences.append(a_support)

            b_support, start_offset = compute_supported_node_sequence(
                start_offset,
                b_length,
                edge_length,
                row['path_align_start'],
                row['path_align_end']
            )
            b_support['name'] = node_b
            b_support['orientation'] = b_orient
            b_support['read_name'] = read_name
            b_support['quality'] = row['mapq']
            b_support['length'] = b_length
            b_support['align_type'] = 'path_member'
            supported_sequences.append(b_support)

    return path_nodes, nodes_per_block, supported_edges, supported_sequences


def compute_pairwise_edges(read_name, node_lengths, path_nodes, supported_edges, alignments):

    if alignments['read_aln_block'].nunique() == 1:
        return supported_edges
    pairs_done = set()
    for block_a, align_a in alignments.groupby('read_aln_block'):
        for block_b, align_b in alignments.groupby('read_aln_block'):
            if block_a == block_b:
                continue
            elif (block_a, block_b) in pairs_done or (block_b, block_a) in pairs_done:
                continue
            else:
                node_iter = itt.product(align_a['path'], align_b['path'])
                mapq_iter = itt.product(align_a['mapq'], align_b['mapq'])
                for (node_a, node_b), (mapq_a, mapq_b) in zip(node_iter, mapq_iter):
                    a_name = node_a[1:]
                    b_name = node_b[1:]
                    if a_name == b_name:
                        continue
                    # if both nodes are members of a path, assume
                    # that is the more informative alignment
                    if a_name in path_nodes and b_name in path_nodes:
                        continue
                    a_orient = ORIENTATION_MAP[node_a[0]]
                    b_orient = ORIENTATION_MAP[node_b[0]]
                    attributes = {
                        'a_length': node_lengths[a_name],
                        'a_orient': a_orient,
                        'b_length': node_lengths[b_name],
                        'b_orient': b_orient,
                        'overlap': 0,
                        'source': 'alignment',
                        'quality': (mapq_a + mapq_b) // 2,
                        'read_name': read_name,
                    }
                    supported_edges.append((a_name, b_name, attributes))
        pairs_done.add((block_a, block_b))
        pairs_done.add((block_b, block_a))
    return supported_edges


def process_singleton_alignments(
    read_name,
    node_lengths,
    path_nodes,
    nodes_per_block,
    supported_edges,
    supported_sequences,
    alignments
):
    """
    Consider cases:
    - true singleton alignment: no overlap, no path member
    - assembly error: no overlap, but path member
    --- does not induce edge as node is already member of a path
    - overlap path alignment, no path member
    --- induces new edge to all path nodes (alternative route)
    - overlaps path alignment and is member of path
    --- count as multi-mapper
    """
    singleton_block_counts = alignments['read_aln_block'].value_counts()
    for _, row in alignments.iterrows():
        node_name = row['path']
        node_orientation = ORIENTATION_MAP[node_name[0]]
        node_name = node_name[1:]
        node_length = node_lengths[node_name]
        node_support, _ = compute_supported_node_sequence(
            0, node_length, 0,
            row['path_align_start'],
            row['path_align_end']
        )
        node_support['name'] = node_name
        node_support['orientation'] = node_orientation
        node_support['read_name'] = read_name
        node_support['quality'] = row['mapq']
        node_support['length'] = node_length
        if row['read_aln_block'] not in nodes_per_block and node_name not in path_nodes:
            # true singleton
            if singleton_block_counts[row['read_aln_block']] == 1:
                node_support['align_type'] = 'singleton'
            else:
                node_support['align_type'] = 'multiple'
        elif row['read_aln_block'] in nodes_per_block and node_name in path_nodes:
            node_support['align_type'] = 'multiple'
        elif row['read_aln_block'] not in nodes_per_block and node_name in path_nodes:
            # probably assembly error, no new edge necessary
            node_support['align_type'] = 'path_disconnect'
        elif row['read_aln_block'] in nodes_per_block and node_name not in path_nodes:
            # bubble / alternate path, induces new edge
            node_support['align_type'] = 'connect'
            for node_b, mapq_b in nodes_per_block[row['read_aln_block']]:
                attributes = {
                    'a_length': node_length,
                    'a_orient': node_orientation,
                    'b_length': node_lengths[node_b],
                    'b_orient': 1,
                    'overlap': 0,
                    'source': 'alignment',
                    'quality': (row['mapq'] + mapq_b) // 2,
                    'read_name': read_name
                }
                supported_edges.append((node_name, node_b, dict(attributes)))
                attributes['b_orient'] = -1
                supported_edges.append((node_name, node_b, dict(attributes)))
        else:
            raise RuntimeError(f'Unexpected program state: {read_name} / {alignments}')
        supported_sequences.append(node_support)

    # connect all singleton alignments of different blocks
    supported_edges = compute_pairwise_edges(
        read_name,
        node_lengths,
        path_nodes,
        supported_edges,
        alignments
    )
    return supported_edges, supported_sequences


def process_read_alignments(node_lengths, edge_lengths, segment_re, args):
    """
    Umbrella function to process all alignments of one
    single read (to an assembly graph) in one go
    - transforms alignments to a path |P| > 1 to a series
    of node-pair alignments taking edge lengths from the
    cached graph data (i.e., edge lengths are not accurate here)
    - attempts to classify alignments that may indicate
    problems in the graph (misassembly, disconnect from path etc.)
    """
    read_name, alignments = args
    alignments['path_size'] = alignments['path'].str.count('>') + alignments['path'].str.count('<')

    # enumerate read alignments by query positions (= alignment coordinates on read)
    # to check for overlapping alignments that may or may not be meaningful
    alignments = enumerate_overlapping_query_blocks(alignments)

    # process all alignments to paths |P| > 1
    aln_to_path = alignments.loc[alignments['path_size'] > 1, :]
    path_nodes, nodes_per_block, supported_edges, supported_sequences = process_path_alignments(
        read_name,
        node_lengths,
        edge_lengths,
        segment_re,
        aln_to_path
    )

    # process all alignments to paths |P| < 2
    aln_to_path = alignments.loc[alignments['path_size'] < 2, :]
    supported_edges, supported_sequences = process_singleton_alignments(
        read_name,
        node_lengths,
        path_nodes,
        nodes_per_block,
        supported_edges,
        supported_sequences,
        aln_to_path
    )

    return supported_edges, supported_sequences


def check_required_cache_files(cache_folder, cache_prefix, required_suffix, required_names):

    if not cache_folder.is_dir():
        logger.error(f'Cache folder does not exist: {cache_folder}')

    err_msg = ''
    required_caches = dict()
    for suffix, name in zip(required_suffix, required_names):
        cache_file = list(cache_folder.glob(f'{cache_prefix}*{suffix}'))
        if not cache_file:
            err_msg += f'Required cache "{name}" missing from folder {cache_folder}\n'
        elif len(cache_file) > 1:
            err_msg += f'Multiple cache files for "{name}": {sorted(cache_file)}\n'
        else:
            required_caches[name] = cache_file[0]
    if err_msg:
        logger.error(err_msg)
        raise RuntimeError('Required graph cache incomplete')
    return required_caches


def prepare_cache_files(cache_folder, cache_prefix, ignore_cache):

    required_cache_suffix = ['.node-lengths.cache.pck', '.edge-lengths.cache.pck']
    required_cache_names = ['node-lengths', 'edge-lengths']
    graph_cache = check_required_cache_files(cache_folder, cache_prefix, required_cache_suffix, required_cache_names)

    file_suffix = [
        '.gaf.cache.h5',
        '.graph-align.cache.pck',
        '.edges-align.cache.pck',
    ]
    cache_names = [
        'gaf-input',
        'graph-align',
        'edges-align',
    ]

    cache_files = dict()
    for s, n in zip(file_suffix, cache_names):
        file_path = pl.Path(cache_folder, cache_prefix + f'{s}')
        cache_files[n] = file_path

    incomplete = True
    cache_missing = [0] * len(cache_names)
    if not ignore_cache:
        for pos, cache_name in enumerate(cache_names):
            cache_path = cache_files[cache_name]
            if not cache_path.is_file():
                cache_missing[pos] = 1
        incomplete = sum(cache_missing) > 0
        if incomplete:
            logger.warning('"ignore-cache" is not set, but cache is incomplete...')
            for is_missing, cache_name in zip(cache_missing, cache_names):
                if is_missing == 0:
                    continue
                logger.warning(f'Missing cache "{cache_name}", expected file {cache_files[cache_name]}')

    return cache_files, incomplete, graph_cache


def get_tag_info(tag):

    def _conv_int(x):
        return int(x.split(':')[-1])

    def _conv_float(x):
        return float(x.split(':')[-1])

    if ':i:' in tag:
        converter = _conv_int
    elif ':f:' in tag:
        converter = _conv_float
    else:
        raise ValueError(f'Unknown data type of tag: {tag}')
    tag_column_name = {
        'NM': 'edit_distance',
        'AS': 'alignment_score',
        'dv': 'seq_divergence',
        'id': 'seq_identity'
    }
    return tag_column_name[tag[:2]], converter


def prepare_gaf_input(gaf_file_input, gaf_file_cache):

    with xopen.xopen(gaf_file_input, 'rt') as gaf:
        while 1:
            inspect_line = gaf.readline().strip()
            if not inspect_line or inspect_line.startswith('#'):
                continue
            break
    columns = inspect_line.split('\t')
    if len(columns) < 12:
        raise ValueError(f'Less than 12 columns in GAF input - invalid format for file {gaf_file_input}')
    # check mapq column
    try:
        mapq = int(columns[11])
    except ValueError:
        raise ValueError(f'MAPQ value (column 12) appears to be non-numeric: {columns[11]}')
    if mapq == 255:
        raise ValueError('MAPQ value in first line is set to 255, i.e., missing. This is not supported.')

    column_names = list(GAF_DEFAULT_COLUMNS)
    use_names = list(GAF_DEFAULT_COLUMNS)
    column_types = GAF_DEFAULT_COLUMN_TYPE_MAP

    def conv_orient(x):
        if x != '+':
            raise ValueError('Illegal read alignment orientation: "+"')
        return 1

    converters = {
        'read_align_orientation': conv_orient
    }
    skip_count = 0
    for c in columns[12:]:
        if c[:2] not in ['NM', 'AS', 'dv', 'id']:
            logger.warning(f'Skipping over unknown tag: {c}')
            column_names.append(f'skip_{skip_count}')
            skip_count += 1
            continue
        col_name, col_converter = get_tag_info(c)
        column_names.append(col_name)
        use_names.append(col_name)
        converters[col_name] = col_converter
    logger.debug(f'Found {len(column_names)} columns in GAF input of which {len(use_names)} will be used.')

    alignments = pd.read_csv(
        gaf_file_input,
        sep='\t',
        header=None,
        index_col=None,
        names=column_names,
        usecols=use_names,
        dtype=column_types,
        converters=converters
    )
    logger.debug('Dumping typed GAF alignments to cache file...')
    with pd.HDFStore(gaf_file_cache, 'w', complevel=9) as hdf:
        hdf.put('cache', alignments)
    return alignments


def read_gaf_input(gaf_file_path, gaf_file_cache, ignore_cache):

    if not ignore_cache and gaf_file_cache.is_file():
        logger.debug(f'Loading cached alignment data from {gaf_file_cache}')
        alignments = pd.read_hdf(gaf_file_cache, 'cache')
    else:
        alignments = prepare_gaf_input(gaf_file_path, gaf_file_cache)
    logger.debug(f'Loaded {alignments.shape[0]} alignment records...')

    read_counter = 0
    for read_name, read_alignments in alignments.groupby('read_name'):
        read_counter += 1
        yield read_name, read_alignments
    logger.debug(f'Iterated through all alignments grouped by {read_counter} reads')
    return


def process_gaf_input(gaf_file_path, gaf_file_cache, ignore_cache, node_lengths, edge_lengths, num_jobs):
    """
    The current implementation only works if the GAF file has a informative MAPQ column
    """
    segment_re = re.compile(GRAPH_SEGMENT_REGEXP)
    process_read_aln = fnt.partial(process_read_alignments, node_lengths, edge_lengths, segment_re)

    processed_records = col.Counter()
    edges = []
    alignments = []
    with mp.Pool(num_jobs) as pool:
        resit = pool.imap_unordered(process_read_aln, read_gaf_input(gaf_file_path, gaf_file_cache, ignore_cache))
        for supported_edges, read_alignments in resit:
            processed_records['reads'] += 1
            if processed_records['reads'] % 10000 == 0:
                logger.info(f'Processed {processed_records["reads"]} GAF records')
            processed_records['edges'] += len(supported_edges)
            processed_records['sequences'] += len(read_alignments)
            edges.extend(supported_edges)
            alignments.extend(read_alignments)
    logger.debug('All GAF records processed')
    logger.debug(f'Total reads {processed_records["reads"]} '
                 f'/ edges {processed_records["edges"]} '
                 f'/ sequence support {processed_records["sequences"]}')

    return alignments, edges


def load_cached_graph_data(graph_cache):

    with open(graph_cache['node-lengths'], 'rb') as cache:
        node_lengths = pck.load(cache)
    logger.debug(f'Loaded cached length data for {len(node_lengths)} graph nodes (segments)')

    with open(graph_cache['edge-lengths'], 'rb') as cache:
        edge_lengths = pck.load(cache)
    logger.debug(f'Loaded cached edge length data for {len(edge_lengths)} graph edges (links)')

    return node_lengths, edge_lengths


def cache_alignment_data(alignments, edges, cache_files):

    logger.debug('Dumping alignment-derived edge / link cache...')
    with open(cache_files['edges-align'], 'wb') as cache:
        pck.dump(edges, cache)

    logger.debug('Dumping processed alignments cache...')
    with open(cache_files['graph-align'], 'wb') as cache:
        pck.dump(alignments, cache)

    return


def run_process_gaf(args):

    logger.debug(f'Processing GAF input file {args.align}')

    logger.debug(f'Preparing cache files in folder {args.cache_folder}')
    cache_files, cache_incomplete, graph_cache = prepare_cache_files(
        args.cache_folder,
        args.cache_prefix,
        args.ignore_cache
    )

    if cache_incomplete or args.ignore_cache:

        logger.debug('Loading cached graph information')
        node_lengths, edge_lengths = load_cached_graph_data(graph_cache)

        alignments, edges = process_gaf_input(
            args.align,
            cache_files['gaf-input'],
            args.ignore_cache,
            node_lengths,
            edge_lengths,
            args.num_jobs
        )
        cache_alignment_data(alignments, edges, cache_files)
    else:
        logger.warning(
            'Graph alignment cache seems complete and "ignore cache" is not set. '
            'There is nothing to do - exiting...'
        )

    return


if __name__ == '__main__':
    parser = process_gaf_cli_parser(None)
    args = parser.parse_args()
    setattr(args, 'num_jobs', 2)
    run_process_gaf(args)
