import collections as col
import argparse as argp
import pathlib as pl
import functools as fnt
import logging as logging
import multiprocessing as mp
import re
import shelve
import pickle as pck

import pandas as pd
import xopen
import networkx as nx

from ontourage import ORIENTATION_MAP
from ontourage.structs import Node, Edge

logger = logging.getLogger()


def process_gfa_cli_parser(module_parsers):

    name = 'process-graph'
    desc = 'Process a GFAv1 input graph and cache its data'

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
        "--graph",
        "--gfa",
        "-g",
        required=True,
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="graph",
        metavar="GFA",
        help="Path to input GFAv1 file. Can be compressed (gzip, bzip, xz)."
    )

    io_group.add_argument(
        "--cache-folder",
        "-cf",
        required=True,
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="cache_folder",
        metavar="CACHE-FOLDER",
        help="Path to a writable folder for caching. Will be created if it does not exist."
    )
    io_group.add_argument(
        "--cache-prefix",
        "-cp",
        default="ONTourage",
        type=str,
        dest="cache_prefix",
        metavar="CACHE-PREFIX",
        help="Specify a prefix for the individual cache files. Default: ONTourage"
    )
    io_group.add_argument(
        "--cache-limit",
        "-l",
        default=536870912,
        type=lambda x: convert_cache_limit(x),
        dest="cache_limit",
        metavar="CACHE-LIMIT",
        help="Specify a (very rough) cache limit for segment sequence caching (in byte, suffix K/M/G recognized). "
             "Default: 512M"
    )
    io_group.add_argument(
        "--ignore-cache",
        "-i",
        action="store_true",
        default=False,
        dest="ignore_cache",
        help="Ignore existing cached data and rerun everything."
    )
    io_group.add_argument(
        "--cache-tables",
        "-ctab",
        action="store_true",
        default=False,
        dest="cache_tables",
        help="In addition to the cache dump in Python's pickle format, dump Edge and Node informaton "
             "also as Pandas DataFrames (HDF files) for easier post-processing with "
             "other tools/libraries. Default: False"
    )

    process_group = parser.add_argument_group('Data processing')

    process_group.add_argument(
        "--no-sequence-stats",
        "-noss",
        default=False,
        action="store_true",
        dest="no_seq_stats",
        help="If the input graph (GFA file) contains node/segment sequences, do not compute "
             "sequence composition statistics (much faster). Default: False"
    )

    parser.set_defaults(execute=run_process_gfa)
    if module_parsers is None:
        return parser
    else:
        return module_parsers


def convert_cache_limit(user_input):

    try:
        cache_limit = int(user_input)
    except ValueError:
        mobj = re.match('[0-9]+', user_input)
        if mobj is None:
            raise ValueError(f'Cannot parse user defined cache limit: {user_input}')
        numeric = int(user_input[:mobj.end()])
        suffix = user_input[mobj.end():].upper()
        factors = {s: 1024 for s in ['K', 'KB', 'KIB']}
        factors.update({s: 1024**2 for s in ['M', 'MB', 'MIB']})
        factors.update({s: 1024**3 for s in ['G', 'GB', 'GIB']})
        try:
            factor = factors[suffix]
        except KeyError:
            raise ValueError(f'Cannot recognize user defined cache limit suffix: {suffix}')
        cache_limit = int(numeric * factor)
    return cache_limit


def get_segment_parser(use_column, no_seq_stats):
    if use_column == 2:
        return fnt.partial(parse_segment_line_has_seq, no_seq_stats)
    else:
        return fnt.partial(parse_segment_line_no_seq, use_column)


def parse_segment_line_has_seq(no_seq_stats, segment_line):

    _, segment_name, segment_sequence = segment_line.split('\t')[:3]
    if no_seq_stats:
        node = Node(segment_name, len(segment_sequence))
    else:
        node = Node(segment_name)
        node.add_sequence_composition(segment_sequence)
    return node, segment_sequence


def parse_segment_line_no_seq(ln_tag_position, segment_line):

    segment_info = segment_line.split('\t')
    segment_name = segment_info[1]
    segment_length = int(segment_info[ln_tag_position].split(':')[-1])
    node = Node(segment_name, segment_length)
    return node, None


def check_gfa_has_sequence(gfa_file_path):

    with xopen.xopen(gfa_file_path, 'rt') as graph:
        while 1:
            inspect_line = graph.readline()
            if not inspect_line.strip() or inspect_line[0] in ['#', 'H']:
                continue
            break
    if inspect_line[0] != 'S':
        raise ValueError(f'Expect first graph entry to be a segment / node / "S", and not: {inspect_line[:5]}')
    inspect_line = inspect_line.strip().split('\t')
    if len(inspect_line) < 3:
        raise ValueError(f'Invalid segment line (must have at least 3 columns): {inspect_line}')

    has_sequence = inspect_line[2] != '*'
    if not has_sequence:
        logger.warning('GFA file does not contain segment sequences - checking for LN tag...')
        # check if we find the LN tag - sequence length is required for certain operations
        tag_ln_column = None
        for pos, item in enumerate(inspect_line):
            if item.startswith('LN:i:'):
                tag_ln_column = pos
                break
        if tag_ln_column is None:
            raise ValueError('No sequence and no LN tag detected - sequence length information is required')
    else:
        logger.debug('GFA file contains sequences')
        tag_ln_column = 2
    return has_sequence, tag_ln_column


def parse_link_line(link_line):

    a_node, a_orient, b_node, b_orient, overlap = link_line.split('\t')[1:6]
    a_orient = ORIENTATION_MAP[a_orient]
    b_orient = ORIENTATION_MAP[b_orient]
    overlap = int(overlap.strip('M'))  # only support "(M)atch" in CIGAR string
    edge = Edge(a_node, a_orient, b_node, b_orient, overlap, -1)
    return edge


def prepare_cache_files(cache_folder, cache_prefix, gfa_has_sequence, ignore_cache):
    cache_folder.mkdir(parents=True, exist_ok=True)

    file_suffix = [
        '.node-lengths.cache.pck',
        '.edge-lengths.cache.pck',
        '.edges-graph.cache.pck',
        '.nodes-graph.cache.pck'
    ]
    cache_names = [
        'node-lengths',
        'edge-lengths',
        'edges-graph',
        'nodes-graph'
    ]
    if gfa_has_sequence:
        file_suffix.append('.sequences.cache')
        cache_names.append('sequences')

    cache_files = dict()
    for s, n in zip(file_suffix, cache_names):
        file_path = pl.Path(cache_folder, cache_prefix + f'{s}')
        cache_files[n] = file_path

    incomplete = True
    cache_missing = [0] * len(cache_names)
    if not ignore_cache:
        for pos, cache_name in enumerate(cache_names):
            cache_path = cache_files[cache_name]
            if cache_name == 'sequences':
                if not gfa_has_sequence:
                    continue
                # shelve module may add different file extensions
                # depending on the underlying db type available
                glob_pattern = f'{cache_prefix}.sequences.cache*'
                any_exist = list(cache_folder.glob(glob_pattern))
                if any_exist:
                    continue
                cache_missing[pos] = 1
            else:
                if not cache_path.is_file():
                    cache_missing[pos] = 1
        incomplete = sum(cache_missing) > 0
        if incomplete:
            logger.warning('"ignore-cache" is not set, but cache is incomplete...')
            for is_missing, cache_name in zip(cache_missing, cache_names):
                if is_missing == 0:
                    continue
                logger.warning(f'Missing cache "{cache_name}", expected file {cache_files[cache_name]}')

    return cache_files, incomplete


def process_gfa_line(seq_column, no_seq_stats, gfa_line):

    if gfa_line[0] == 'S':
        segment_parser = get_segment_parser(seq_column, no_seq_stats)
        node, node_sequence = segment_parser(gfa_line)
        return node, node_sequence
    elif gfa_line[0] == 'L':
        edge = parse_link_line(gfa_line)
        return edge, None
    else:
        raise ValueError(f'Unknown GFA line type {gfa_line.strip()}')


def read_gfa_input(gfa_file_path):

    with xopen.xopen(gfa_file_path, 'rt') as graph:
        for line in graph:
            if line[0] not in ['S', 'L']:
                continue
            yield line
    return


def cache_segment_sequences(cache_file, cache_flag, sequences):

    with shelve.open(str(cache_file), flag=cache_flag) as cache:
        for k, v in sequences.items():
            cache[k] = v
    return


def process_gfa_input(gfa_file_path, seq_column, no_seq_stats, seq_cache_file, num_jobs, cache_limit):
    """
    The current implementation only works if segments precede links
    in the input GFA (to add the segment length info to the links)
    """

    process_line = fnt.partial(process_gfa_line, seq_column, no_seq_stats)

    segment_lengths = dict()
    edge_lengths = dict()
    sequences = dict()
    edges = []
    nodes = []
    cached_sequence_length = 0
    sequence_cache_flag = 'c'
    processed_records = col.Counter()

    with mp.Pool(num_jobs) as pool:
        resit = pool.imap_unordered(process_line, read_gfa_input(gfa_file_path))
        for item, node_seq in resit:
            processed_records['total'] += 1
            if processed_records['total'] % 1000 == 0:
                logger.debug(f'Processed {processed_records["total"]} GFA records')
            if isinstance(item, Node):
                processed_records['segment/node'] += 1
                segment_lengths[item.name] = item.length
                if node_seq is not None:
                    sequences[item.name] = node_seq
                    cached_sequence_length += len(node_seq)
                if cached_sequence_length > cache_limit:
                    logger.debug(f'Caching sequences (~{cached_sequence_length//1024**2} MB)...')
                    cache_segment_sequences(seq_cache_file, sequence_cache_flag, sequences)
                    sequences = dict()
                    cached_sequence_length = 0
                    sequence_cache_flag = 'a'
                nodes.append(item)
            elif isinstance(item, Edge):
                processed_records['link/edge'] += 1
                edges.append(item)
                edge_lengths[item._id] = item.length
            else:
                item_type = type(item)
                raise ValueError(f'Unexpected GFA line type: {item_type}')
    logger.debug(f"GFA processed - {processed_records['segment/node']} segments and {processed_records['link/edge']} links")
    if cached_sequence_length > 0:
        logger.debug(f'Caching remaining sequences (~{cached_sequence_length//1024**2} MB)...')
        cache_segment_sequences(seq_cache_file, sequence_cache_flag, sequences)
    return nodes, edges, segment_lengths, edge_lengths


def cache_graph_data(nodes, node_lengths, edges, edge_lengths, cache_files, dump_tables):

    logger.debug('Dumping node / segment cache...')
    with open(cache_files['nodes-graph'], 'wb') as cache:
        pck.dump(nodes, cache)
    if dump_tables:
        logger.debug('Dumping node / segment cache as HDF/table...')
        hdf_path = cache_files['nodes-graph'].with_suffix('.h5')
        with pd.HDFStore(hdf_path, 'w', complevel=9) as hdf:
            dump_df = pd.DataFrame.from_records([n.as_dict() for n in nodes])
            hdf.put('cache', dump_df, format='fixed')

    logger.debug('Dumping edge / link cache...')
    with open(cache_files['edges-graph'], 'wb') as cache:
        pck.dump(edges, cache)
    if dump_tables:
        logger.debug('Dumping edge / link cache as HDF/table...')
        hdf_path = cache_files['edges-graph'].with_suffix('.h5')
        with pd.HDFStore(hdf_path, 'w', complevel=9) as hdf:
            dump_df = pd.DataFrame.from_records([e.as_dict() for e in edges])
            hdf.put('cache', dump_df, format='fixed')

    logger.debug('Dumping node / segment length cache...')
    with open(cache_files['node-lengths'], 'wb') as cache:
        pck.dump(node_lengths, cache)

    logger.debug('Dumping edge / link length cache...')
    with open(cache_files['edge-lengths'], 'wb') as cache:
        pck.dump(edge_lengths, cache)
    return


def annotate_graph_nodes(nodes, edges):
    """
    Annotate graph nodes with their connected component
    id and cardinality, and their local node betweenness centrality
    (local = within the connected component).
    For this annotation/computation, the edge direction
    is ignored
    """
    logger.debug('Annotating nodes in graph')
    g = nx.Graph()
    logger.debug('Adding edges to graph...')
    g.add_edges_from([(e.node_a, e.node_b) for e in edges])
    logger.debug('Computing node betweeness centrality...')
    cc_members = dict()
    node_centralities = dict()
    cc_sizes = dict()
    for cc_id, cc_nodes in enumerate(nx.connected_components(g), start=0):
        subgraph = g.subgraph(cc_nodes).copy()
        centralities = nx.betweenness_centrality(subgraph, normalized=True, endpoints=True)
        cc_members.update(dict((n, cc_id) for n in cc_nodes))
        cc_sizes[cc_id] = len(cc_nodes)
        node_centralities.update(centralities)
    logger.debug('Computing node degrees...')
    degrees = col.defaultdict(col.Counter)
    for e in edges:
        degrees[e.node_a][('out', e.a_orientation)] += 1
        degrees[e.node_b][('in', e.b_orientation)] += 1
    logger.debug('Updating node information...')
    # networkx does not count isolated nodes as connected components
    # manually add CC IDs for these cases
    manual_cc_id = cc_id + 1
    for n in nodes:
        n.set_node_degree(degrees[n.name])
        try:
            n.cc_id = cc_members[n.name]
        except KeyError:
            n.cc_id = manual_cc_id
            n.cc_cardinality = 1
            n.centrality = 1
            manual_cc_id += 1
        else:
            n.cc_cardinality = cc_sizes[n.cc_id]
            n.centrality = node_centralities[n.name]
    logger.debug('Annotation of graph nodes completed')
    return


def run_process_gfa(args):

    logger.debug(f'Processing GFA input file {args.graph}')
    logger.debug('Checking if GFA contains segment sequences...')
    gfa_has_sequence, seq_column = check_gfa_has_sequence(args.graph)

    logger.debug(f'Preparing cache files in folder {args.cache_folder}')
    cache_files, cache_incomplete = prepare_cache_files(
        args.cache_folder,
        args.cache_prefix,
        gfa_has_sequence,
        args.ignore_cache
    )

    if not gfa_has_sequence:
        seq_cache = None
    else:
        seq_cache = cache_files['sequences']

    if args.ignore_cache or cache_incomplete:
        nodes, edges, node_lengths, edge_lengths = process_gfa_input(
            args.graph,
            seq_column,
            args.no_seq_stats,
            seq_cache,
            args.num_jobs,
            args.cache_limit
        )
        annotate_graph_nodes(nodes, edges)
        cache_graph_data(nodes, node_lengths, edges, edge_lengths, cache_files, args.cache_tables)
    else:
        logger.warning(
            'Graph cache seems complete and "ignore cache" is not set. '
            'There is nothing to do - exiting...'
        )
    return


if __name__ == '__main__':
    parser = process_gfa_cli_parser(None)
    args = parser.parse_args()
    setattr(args, 'num_jobs', 2)
    run_process_gfa(args)
