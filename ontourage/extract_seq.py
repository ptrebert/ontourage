
import logging
import pathlib as pl
import argparse as argp
import pickle as pkl
import collections as col
import itertools as itt
import shelve
import io

import networkx as nx
import pandas as pd

from ontourage.constants import REVCOMP_TABLE
from ontourage.structs import Edge


logger = logging.getLogger()


def extract_seq_cli_parser(module_parsers):

    name = 'extract-seq'
    desc = 'Extract the genomic sequence between two nodes in a graph.'

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
        "--cache-folder",
        "-cf",
        required=True,
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="cache_folder",
        metavar="CACHE-FOLDER",
        help="Path to an existing folder for caching. "
             "Note that the folder has to be the same as used for the graph (GFA) caching."
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
        "--path-nodes",
        "-pn",
        type=str,
        nargs="+",
        dest="path_nodes",
        metavar="PATH-NODE-SPEC",
        help="Specify sequences to extract as path node names: START,[STEP,[STEP ...]],END "
             "If steps are specified, the path from start to end will be decomposed "
             "to (forcibly) include all steps. In that case, the direct path from start "
             "to end may be shorter even if 'shortest' is specified as path finding strategy. "
             "This argument can also be provided as row-oriented text file, with one line "
             "specifying one path as START,END (plus arbitrarily many steps, if applicable)"
    )

    io_group.add_argument(
        "--output-fasta",
        "-of",
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="output_fasta",
        required=True,
        metavar="OUTPUT-FASTA",
        help="Specify the path to the output FASTA file. If the option '--split-fasta' is set, this "
             "is interpreted as a folder path."
    )

    io_group.add_argument(
        "--split-fasta",
        "-sf",
        action="store_true",
        default=False,
        dest="split_fasta",
        help="Write one FASTA per extracted path sequence. In that case, '--output-fasta' is "
             "interpreted as an output folder. Default: False"
    )

    process_group = parser.add_argument_group('Data processing')

    process_group.add_argument(
        "--path-finding",
        "-pf",
        default="shortest",
        choices=["shortest"],
        type=str,
        dest="path_finding",
        help="Specify strategy to determine the path between start and end node. "
             "Choices: shortest "
             "Default: shortest"
    )

    process_group.add_argument(
        "--connect-sinks",
        "-cs",
        action="store_true",
        default=False,
        dest="connect_sinks",
        help="If a node is a sink (out degree zero, in degree non-zero), invert all "
             "incoming edges and add these to the graph, too. Default: False"
    )

    process_group.add_argument(
        "--fix-path-direction",
        "-fpd",
        action="store_true",
        default=False,
        dest="fix_path_direction",
        help="If set, consider directionality of the path direction as fixed/strict, i.e. "
             "if there is no path from A to B, do not check for a path from B to A. "
             "Has no effect when '--undirected-graph' is set. "
             "Default: False"
    )

    process_group.add_argument(
        "--undirected-graph",
        "-udg",
        action="store_true",
        default=False,
        dest="undirected_graph",
        help="Build and traverse the graph ignoring edge directionality. "
             "Strand sense of the node sequences will still be considered. "
             "Default: False"
    )

    process_group.add_argument(
        "--allow-skipping-steps",
        "-skip",
        action="store_true",
        default=False,
        dest="skip_steps",
        help="When traversing the graph along a path with dead-ends, allow to skip "
             "to the next step instead of aborting the traversal. "
             "Default: False"
    )

    parser.set_defaults(execute=run_extract_seq)
    if module_parsers is None:
        return parser
    else:
        return module_parsers


def parse_path_specification_entry(path_spec):

    path_components = path_spec.strip().split(',')
    if len(path_components) < 2:
        raise ValueError(f"Valid path specification must have at least start and end node: {path_spec}")
    start = path_components[0]
    end = path_components[-1]
    logger.debug(f"Extracting path/sequence: {' -> '.join(path_components)}")
    return (start, end), path_components


def iter_path_spec_file(file_path):

    with open(file_path, 'r') as listing:
        for line in listing:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            seq_name, extract_seq = parse_path_specification_entry(line)
            yield seq_name, extract_seq
    return


def iter_path_spec_listing(param_listing):

    for record in param_listing:
        seq_name, extract_seq = parse_path_specification_entry(record)
        yield seq_name, extract_seq
    return


def get_path_spec_iterator(arg_path_spec):

    try:
        logger.debug('Checking if path spec is given as file...')
        path_spec = pl.Path(arg_path_spec[0]).resolve().absolute()
        if not path_spec.is_file():
            raise TypeError(f"This is not an existing file: {arg_path_spec[0]}")
        logger.debug('Checking if path spec file contains at least one valid entry...')
        for _ in iter_path_spec_file(path_spec):
            break
        path_spec_iter = iter_path_spec_file
    except TypeError:
        path_spec = arg_path_spec
        logger.debug('Assuming path spec is provided as command-line listing...')
        logger.debug('Checking for at least one valid entry...')
        for _ in iter_path_spec_listing(arg_path_spec):
            break
        path_spec_iter = iter_path_spec_listing
    return path_spec_iter, path_spec


def check_required_cache_files(cache_folder, cache_prefix, required_suffix, required_names):

    if not cache_folder.is_dir():
        logger.error(f'Cache folder does not exist: {cache_folder}')

    err_msg = ''
    required_caches = dict()
    for suffix, name in zip(required_suffix, required_names):
        if name == 'sequences':
            cache_file = list(cache_folder.glob(f'{cache_prefix}*{suffix}*'))
        else:
            cache_file = list(cache_folder.glob(f'{cache_prefix}*{suffix}'))
        if not cache_file:
            err_msg += f'Required cache "{name}" missing from folder {cache_folder}\n'
        elif len(cache_file) > 1 and name != 'sequences':
            err_msg += f'Multiple cache files for "{name}": {sorted(cache_file)}\n'
        else:
            if name == 'sequences':
                required_caches[name] = cache_folder / pl.Path(f'{cache_prefix}{suffix}')
            else:
                required_caches[name] = cache_file[0]
    if err_msg:
        logger.error(err_msg)
        raise RuntimeError('Required graph cache incomplete')
    return required_caches


def is_sink_node(node):
    return node.deg_total > 0 and node.deg_out_plus == 0 and node.deg_out_minus == 0


def build_directed_graph(graph_node_cache, graph_edge_cache, connect_sinks):

    graph = nx.DiGraph()
    with open(graph_edge_cache, 'rb') as dump:
        graph_edges = pkl.load(dump)
    logger.debug(f"Loaded {len(graph_edges)} edges from cache")
    if not connect_sinks:
        graph.add_edges_from([e.as_nx_edge() for e in graph_edges])
    else:
        logger.debug("Fully connecting sink nodes...")
        logger.debug("Loading cached node data...")
        with open(graph_node_cache, 'rb') as dump:
            graph_nodes = {n.name: n for n in pkl.load(dump)}
        for edge in graph_edges:
            a, b, attr = edge.as_nx_edge()
            if is_sink_node(graph_nodes[edge.node_b]):
                inv_a, inv_b, inv_attr = edge.flip().as_nx_edge()
                logger.debug(f"Connecting sink via flip: {inv_a} => {inv_b} | {inv_attr}")
                graph.add_edge(inv_a, inv_b, **inv_attr)
                if inv_a[-1] != inv_b[-1]:
                    cross_a, cross_b, cross_attr = edge.cross_flip().as_nx_edge()
                    logger.debug(f"Connecting sink via cross-flip: {cross_a} => {cross_b} | {cross_attr}")
                    graph.add_edge(cross_a, cross_b, **cross_attr)
            graph.add_edge(a, b, **attr)
    logger.debug(f'Directected Graph built complete. Graph size: {graph.size()}')
    return graph


def needs_self_connect(node):
    """
    For undirected graphs
    """
    not_singleton = node.deg_total > 0
    plus_in_no_out = (node.deg_in_plus > 0) & (node.deg_out_plus == 0)
    plus_out_no_in = (node.deg_in_plus == 0) & (node.deg_out_plus > 0)
    minus_in_no_out = (node.deg_in_minus > 0) & (node.deg_out_minus == 0)
    minus_out_no_in = (node.deg_in_minus == 0) & (node.deg_out_minus > 0)
    blocked_source = plus_out_no_in & minus_out_no_in
    blocked_sink = plus_in_no_out & minus_in_no_out
    return not_singleton & (blocked_source | blocked_sink)


def build_undirected_graph(graph_node_cache, graph_edge_cache, connect_sinks):

    graph = nx.Graph()
    with open(graph_edge_cache, 'rb') as dump:
        graph_edges = pkl.load(dump)
    if not connect_sinks:
        graph.add_edges_from([e.as_nx_edge() for e in graph_edges])
    else:
        logger.debug('Checking if orientation loops are needed for undirected graph')
        logger.debug('Loading cached node data...')
        with open(graph_node_cache, 'rb') as dump:
            graph_nodes = {n.name: n for n in pkl.load(dump)}
        for edge in graph_edges:
            a, b, attr = edge.as_nx_edge()
            if needs_self_connect(graph_nodes[edge.node_a]):
                logger.debug(f"Adding orientation self-loop for A node {edge.node_a}")
                tmp_edge = Edge(edge.node_a, 1, edge.node_a, -1, 0, 0, edge_type='loop')
                loop_a, loop_b, loop_attr = tmp_edge.as_nx_edge()
                graph.add_edge(loop_a, loop_b, **loop_attr)
            if needs_self_connect(graph_nodes[edge.node_b]):
                logger.debug(f"Adding orientation self-loop for B node {edge.node_b}")
                tmp_edge = Edge(edge.node_b, 1, edge.node_b, -1, 0, 0, edge_type='loop')
                loop_a, loop_b, loop_attr = tmp_edge.as_nx_edge()
                graph.add_edge(loop_a, loop_b, **loop_attr)
            graph.add_edge(a, b, **attr)
    logger.debug(f'Undirected graph built complete. Graph size: {graph.size()}')
    return graph


def node_has_orientation(node_name):
    return node_name[-1] in ['+', '-']


def add_node_orientation(node_a, node_b, fix_path_dir, start_orient):
    if node_has_orientation(node_a) and node_has_orientation(node_b):
        oriented_nodes = [(node_a, node_b)]
    elif node_has_orientation(node_a):
        oriented_nodes = [(node_a, f'{node_b}+'), (node_a, f'{node_b}-')]
    elif node_has_orientation(node_b) and start_orient:
        oriented_nodes = [(f'{node_a}{start_orient}', node_b)]
    elif node_has_orientation(node_b):
        oriented_nodes = [(f'{node_a}+', node_b), (f'{node_a}-', node_b)]
    else:
        if start_orient is None:
            oriented_nodes = [
                (f'{node_a}+', f'{node_b}+'),
                (f'{node_a}-', f'{node_b}-'),
                (f'{node_a}+', f'{node_b}-'),
                (f'{node_a}-', f'{node_b}+')
            ]
        else:
            oriented_nodes = [
                (f'{node_a}{start_orient}', f'{node_b}+'),
                (f'{node_a}{start_orient}', f'{node_b}-'),

            ]
    if not fix_path_dir:
        inv_oriented_nodes = [(b, a) for (a, b) in oriented_nodes]
        oriented_nodes.extend(inv_oriented_nodes)
    return oriented_nodes


def concatenate_node_sequences(
    graph,
    grouped_path_nodes,
    seq_cache_file,
    global_start,
    global_end,
    path_is_broken,
    path_infix,
    global_fasta_buffer,
    global_table_buffer
):
    path_state = 'partial' if path_is_broken else 'complete'
    header_infix = '_' if not path_infix else f'_{path_infix}_'
    global_path_name = f'PATH_{path_state}{header_infix}FROM_{global_start}_TO_{global_end}'
    with shelve.open(str(seq_cache_file), 'r') as cache:

        for node_group in grouped_path_nodes:
            concat_end = 0
            concat_seq = ''
            offset = 0
            norm_nodes = [(n[:-1], n[-1], n) for n in node_group]
            group_start = norm_nodes[0][2]
            group_end = norm_nodes[-1][2]
            for (node1, orient1, no1), (node2, orient2, no2) in nx.utils.pairwise(norm_nodes):
                try:
                    seq1 = cache[node1]
                    if orient1 == '-':
                        seq1 = seq1[::-1].translate(REVCOMP_TABLE)
                except KeyError:
                    logger.error(f'Node {node1} / {orient1} sequence is not in cache file {seq_cache_file}')
                    raise KeyError(f'Unknown node / no sequence: {node1}')
                try:
                    seq2 = cache[node2]
                    if orient2 == '-':
                        seq2 = seq2[::-1].translate(REVCOMP_TABLE)
                except KeyError:
                    logger.error(f'Node {node2} / {orient2} sequence is not in cache file {seq_cache_file}')
                    raise KeyError(f'Unknown node / no sequence: {node2}')
                concat_seq += seq1[offset:]
                concat_end = len(concat_seq)
                global_table_buffer.append((
                    global_path_name, group_start, group_end, no1, len(seq1), offset, concat_end - len(seq1), concat_end
                ))
                edge_data = graph.get_edge_data(no1, no2)
                edge_length = edge_data['length']
                if edge_length < 0:
                    raise ValueError(f'Edge has negative length / distance information: {no1} => {no2} | {edge_data}')
                offset = edge_length

            concat_seq += seq2[offset:]
            concat_end = len(concat_seq)
            global_table_buffer.append((
                global_path_name, group_start, group_end, no2, len(seq2), offset, concat_end - len(seq2), concat_end
            ))
            logger.debug(f'Sequence concat complete for path component {group_start} - {group_end}')
            _ = global_fasta_buffer.write(f'>{global_path_name}_{group_start}_{group_end}\n')
            _ = global_fasta_buffer.write(concat_seq + '\n')

    return


def group_path_sequences(stored_nodes, path_is_broken, num_nodes):

    path_sequence = []
    path_collections = []

    while 1:
        try:
            node = stored_nodes.popleft()
        except IndexError:
            break
        if node is not None:
            path_sequence.append(node)
        else:
            if path_sequence:
                path_collections.append(tuple(path_sequence))
                path_sequence = []
    if path_sequence:
        path_collections.append(tuple(path_sequence))

    nodes_processed = sum(map(len, path_collections))
    if nodes_processed != num_nodes:
        raise ValueError(f'Nodes lost during grouping: input {num_nodes} / processed {nodes_processed}')
    path_collections.append(num_nodes)
    path_collections.append(path_is_broken)
    return path_collections


def extract_single_path_sequence(graph, path_nodes, fix_path_dir, skip_steps, start_orient):

    last_step_complete = False
    path_is_broken = False
    node_sequence = col.deque()
    queued_nodes = 0

    init_orientation = start_orient

    for node_a, node_b in nx.utils.pairwise(path_nodes):
        logger.debug(f'Searching for path connecting pair {node_a} / {node_b}')
        for ao, bo in add_node_orientation(node_a, node_b, fix_path_dir, init_orientation):
            # if we found a path up to node_a in the last iteration,
            # make sure to start the traversal at the identical node now
            # for proper continuity of the path sequence
            if last_step_complete:
                last_node_a = node_sequence.pop()
                if last_node_a != ao:
                    node_sequence.append(last_node_a)
                    continue
                node_sequence.append(last_node_a)
            found_component = False
            try:
                # TODO allow other strategies
                logger.debug(f'Searching path: {ao} => {bo}')
                path = nx.shortest_path(graph, ao, bo)
            except nx.NetworkXNoPath:
                logger.warning(f'No path found between nodes {ao} and {bo}')
                continue
            except nx.NodeNotFound:
                logger.warning(f'(Oriented) node not present in graph: {ao} or {bo}')
                continue
            else:
                found_component = True
                logger.debug(f'Found path connecting {ao} and {bo} of length {len(path)}')
                if not last_step_complete:
                    logger.debug('Last step was not completed')
                    node_sequence.append(None)
                    [node_sequence.append(p) for p in path]
                    queued_nodes += len(path)
                else:
                    # NB: it is checked above that we start from
                    # the previous last node
                    [node_sequence.append(p) for p in path[1:]]
                    queued_nodes += len(path) - 1
                break  # break inner loop

        if not found_component:
            last_step_complete = False
            path_is_broken = True
            logger.warning(f'Could not find a path between nodes {node_a} and {node_b}')
            if not skip_steps:
                logger.warning('Skipping steps is not allowed, aborting traversal...')
                logger.debug(f'Emptying node sequence queue of size {queued_nodes}')
                node_sequence.clear()
                break  # break outer loop
        else:
            last_step_complete = True
            init_orientation = bo[-1]

    if node_sequence:
        logger.debug('Grouping path components')
        path_collection = group_path_sequences(node_sequence, path_is_broken, queued_nodes)
    else:
        raise nx.NetworkXNoPath(f'No path found - last pair: {ao} - {bo}')

    return path_collection


def determine_better_path(path_plus, path_minus):
    if path_plus is None and path_minus is None:
        logger.debug('Both paths are empty')
        selected = [None, None]
    elif path_plus is not None and path_minus is None:
        logger.debug('Selecting forward path as final / is not None')
        selected = [path_plus, None]
    elif path_plus is None and path_minus is not None:
        logger.debug('Selecting reverse path as final / is not None')
        selected = [None, path_minus]
    elif path_plus[-1] and not path_minus[-1]:
        # minus is not broken
        logger.debug('Selecting reverse path as final / is complete')
        selected = [None, path_minus]
    elif path_minus[-1] and not path_plus[-1]:
        # plus is not broken
        logger.debug('Selecting forward path as final / is complete')
        selected = [path_plus, None]
    elif path_plus[-2] > path_minus[-2]:
        logger.debug('Selecting forward path as final / is larger')
        selected = [path_plus, None]
    elif path_plus[-2] < path_minus[-2]:
        logger.debug('Selecting reverse path as final / is larger')
        selected = [None, path_minus]
    else:
        # both paths have the same size, check if they cover different
        # node sets
        # strip off orientation symbol from node name
        nodes_frw = sorted(set([n[:-1] for n in itt.chain.from_iterable(path_plus[:-2])]))
        nodes_rev = sorted(set([n[:-1] for n in itt.chain.from_iterable(path_minus[:-2])]))
        if nodes_frw == nodes_rev:
            logger.debug('Both paths identical - selecting forward')
            selected = [path_plus, None]
        else:
            logger.debug('Paths cover different nodes - selecting both')
            selected = [path_plus, path_minus]
    return selected


def extract_path_sequences(cache_files, path_spec_iter, path_spec, seq_cache_file, undirected_graph, connect_sinks, fix_path_dir, skip_steps):

    logger.debug("Building graph from cache...")
    if undirected_graph:
        graph = build_undirected_graph(cache_files['nodes-graph'], cache_files['edges-graph'], connect_sinks)
    else:
        graph = build_directed_graph(cache_files['nodes-graph'], cache_files['edges-graph'], connect_sinks)

    global_fasta_buffer = io.StringIO()
    global_table_buffer = []

    dir_symbol = " ==>> " if fix_path_dir else " <<==>> "
    for (start, end), path_nodes in path_spec_iter(path_spec):
        logger.debug(f"Attempting traversal for path {start}{dir_symbol}{end}")
        if node_has_orientation(start):
            try:
                complete_path = extract_single_path_sequence(
                    graph, path_nodes,
                    fix_path_dir, skip_steps, None
                )
                traversals = [complete_path, None]
            except nx.NetworkXNoPath:
                logger.warning(f"No path found starting from node {start}")
                traversals = [None, None]
        else:
            try:
                complete_path_plus = extract_single_path_sequence(
                    graph, path_nodes,
                    fix_path_dir, skip_steps, '+'
                )
            except nx.NetworkXNoPath:
                logger.warning(f"No path found starting from node {start} in forward orientation")
                complete_path_plus = None
            try:
                complete_path_minus = extract_single_path_sequence(
                    graph, path_nodes,
                    fix_path_dir, skip_steps, '-'
                )
            except nx.NetworkXNoPath:
                logger.warning(f"No path found starting from node {start} in reverse orientation")
                complete_path_minus = None
            traversals = determine_better_path(complete_path_plus, complete_path_minus)

        if any(t is not None for t in traversals):
            logger.debug(f"Generating sequence for path {start}{dir_symbol}{end}")
            if traversals[0] is not None:
                concatenate_node_sequences(
                    graph, traversals[0][:-2], seq_cache_file, start, end, traversals[0][2], 1,
                    global_fasta_buffer, global_table_buffer
                )

            if traversals[1] is not None:
                concatenate_node_sequences(
                    graph, traversals[1][:-2], seq_cache_file, start, end, traversals[1][2], 2,
                    global_fasta_buffer, global_table_buffer
                )

    return global_fasta_buffer, global_table_buffer


def run_extract_seq(args):

    logger.debug("Checking if graph cache seems to be complete")
    # we need edge and sequence information
    required_suffixes = [
        '.edges-graph.cache.pck',
        '.nodes-graph.cache.pck',
        '.sequences.cache'
    ]
    required_names = ['edges-graph', 'nodes-graph', 'sequences']
    cache_files = check_required_cache_files(
        args.cache_folder, args.cache_prefix,
        required_suffixes, required_names
    )
    logger.debug("Checking path spec information...")
    path_spec_iter, path_spec = get_path_spec_iterator(args.path_nodes)

    if args.undirected_graph:
        logger.debug('Forcing "fix path direction" to True because "undirected graph" is set.')
        setattr(args, 'fix_path_direction', True)

    fasta_buffer, table_buffer = extract_path_sequences(
        cache_files,
        path_spec_iter,
        path_spec,
        cache_files['sequences'],
        args.undirected_graph,
        args.connect_sinks,
        args.fix_path_direction,
        args.skip_steps
    )

    with open(args.output_fasta, 'w') as dump:
        _ = dump.write(fasta_buffer.getvalue())
    logger.debug('FASTA dumped')

    path_table = pd.DataFrame.from_records(
        table_buffer,
        columns=['path', 'node', 'group_start', 'group_end', 'seq_length', 'offset', 'seq_start', 'seq_end']
    )
    logger.debug(f'Created path table of size: {path_table.shape}')
    path_table.to_csv(args.output_fasta.with_suffix('.paths.tsv'), sep='\t', header=True, index=False)

    return


if __name__ == '__main__':
    parser = extract_seq_cli_parser(None)
    args = parser.parse_args()
    setattr(args, 'num_jobs', 2)
    run_extract_seq(args)
