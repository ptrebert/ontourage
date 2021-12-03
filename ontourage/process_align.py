import collections as col
import argparse as argp
import pathlib as pl
import functools as fnt
import logging as logging
import itertools as itt
import multiprocessing as mp
import re
import pickle as pck
import hashlib

import pandas as pd
import xopen

from ontourage import GAF_DEFAULT_COLUMN_TYPE_MAP, ORIENTATION_MAP, GAF_DEFAULT_COLUMNS, SEGMENT_ORIENTATION_SYMBOLS
from ontourage.structs import Edge, NodeSupport, IssueRecord


logger = logging.getLogger()


def process_gaf_cli_parser(module_parsers):

    name = 'process-align'
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
    io_group.add_argument(
        "--cache-tables",
        "-ctab",
        action="store_true",
        default=False,
        dest="cache_tables",
        help="In addition to the cache dump in Python's pickle format, dump Edge and NodeSupport informaton "
             "also as Pandas DataFrames (HDF files) for easier post-processing with "
             "other tools/libraries. Default: False"
    )
    io_group.add_argument(
        "--stat-summary",
        "-sts",
        default=None,
        type=str,
        dest="stat_summary",
        metavar="STATISTICS",
        help="Specify a path to save the summary statistic output (key-value pairs). If omitted, "
             "the file will be created one level above the 'cache folder' using the 'cache prefix' "
             "plus the extension '.stats.tsv' as file name: Default: <empty>"
    )
    io_group.add_argument(
        "--issue-records",
        "-irs",
        default=None,
        type=str,
        dest="issue_records",
        metavar="ISSUES",
        help="Specify a path to save the gzipped TSV table of potential issues revealed by the alignments. "
             "If omitted, the file will be created one level above the 'cache folder' using the "
             "'cache prefix' plus the extension '.issues.tsv.gz' as file name: Default: <empty>"
    )


    parser.set_defaults(execute=run_process_gaf)
    if module_parsers is None:
        return parser
    else:
        return module_parsers


def collect_node_sequence_support(read_name, node_lengths, edge_lengths, path_members, aln_record):
    """
    Iteratively collect node sequence support for all path members

    """
    nodes_support = []
    edges_support = []

    quality = aln_record['mapq']
    align_start_in_path = aln_record['path_align_start']
    align_end_in_path = aln_record['path_align_end']
    align_id = aln_record['_id']

    if len(path_members) == 1:
        node_name, node_orient = path_members[0]
        # obviously, no edge support can exist
        node_length = node_lengths[node_name]
        node_support, _ = compute_supported_node_sequence(
            node_name, node_orient,
            read_name, quality,
            0, node_length, 0,
            align_start_in_path, align_end_in_path,
            align_id
        )
        nodes_support.append(node_support)
    else:
        path_offset = 0
        # more elegant in Python 3.10+: itt.pairwise(path_members)
        for (node_a, a_orient), (node_b, b_orient) in zip(path_members[:-1], path_members[1:]):
            edge = Edge(
                node_a, a_orient,
                node_b, b_orient,
                0, quality, read_name,
                align_id
            )
            edge_length = edge_lengths[edge._id]
            edge.length = edge_length

            # note here: for the path U-V-W, edge length 0 is used
            # to compute the new offset when handling node V,
            # because the the correct offset can only be computed
            # using the Edge(V,W), which is not yet available
            # when considering the pair/edge (U,V) in the path.
            if path_offset != 0:
                path_offset -= edge_length
            if path_offset == 0:
                # avoid storing support twice
                a_length = node_lengths[node_a]
                a_support, new_offset = compute_supported_node_sequence(
                    node_a, a_orient, read_name, quality,
                    path_offset, a_length, edge_length,
                    align_start_in_path, align_end_in_path, align_id
                )
                path_offset = new_offset
                nodes_support.append(a_support)

            b_length = node_lengths[node_b]
            b_support, new_offset = compute_supported_node_sequence(
                node_b, b_orient, read_name, quality,
                path_offset, b_length, 0,  # edge length zero - see explanation above
                align_start_in_path, align_end_in_path, align_id
            )
            path_offset = new_offset
            nodes_support.append(b_support)
            edges_support.append(edge)
    return nodes_support, edges_support


def compute_supported_node_sequence(
    node_name,
    node_orient,
    read_name,
    quality,
    start_offset,
    node_length,
    edge_length,
    aln_start_in_path,
    aln_end_in_path,
    aln_id
):
    """
    Compute sequence position of node in path that are supported by a read alignment

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
    assert node_support_end > node_support_start, f"{node_name} - {read_name} - {aln_id} - {node_support_start} - {node_support_end} - {node_start_in_path}"

    support_length = node_support_end - node_support_start
    support_pct = round(support_length / node_length * 100, 2)

    new_offset = node_end_in_path - edge_length

    node_support = NodeSupport(
        node_name,
        node_orient,
        node_length,
        read_name,
        quality,
        node_support_start,
        node_support_end,
        support_length,
        support_pct,
        aln_id
    )

    return node_support, new_offset


def extract_path_members(orient_split, path):
    members = itt.filterfalse(lambda x: not x, orient_split.split(path))
    norm_members = []
    while 1:
        try:
            node_orient = ORIENTATION_MAP[next(members)]
            node_name = next(members)
            norm_members.append((node_name, node_orient))
        except (StopIteration, RuntimeError):
            break
    return norm_members


def check_intra_node_overlap(node_a, node_b):
    """
    check if alignments of the same read hit overlapping
    parts of the same node (or node_a and _b are different)
    """
    same_node = True
    if node_a.name != node_b.name:
        same_node = False
        same_orientation = False
        overlap_bp = 0
        overlap_pct = 0
    else:
        overlap_bp = min(node_a.support_end, node_b.support_end) - max(node_a.support_start, node_b.support_start)
        # if overlap is negative (= alignments are apart), percent overlap should be 0
        overlap_pct = max(0, round(overlap_bp / node_a.length * 100, 2))
        same_orientation = node_a.orientation == node_b.orientation
    return same_node, same_orientation, overlap_bp, overlap_pct


def handle_self_alignments(
    last_node, last_rstart, last_rend,
    this_node, this_rstart, this_rend,
    orient_match,
    ovl_bp_rspace, ovl_pct_rspace,
    ovl_bp_nspace, ovl_pct_nspace,
):
    issue_record = None
    stats_record = None
    put_back = False

    record_last = last_node.as_prefixed_record('a', last_rstart, last_rend)
    record_first = this_node.as_prefixed_record('b', this_rstart, this_rend)
    overlaps = [ovl_bp_rspace, ovl_pct_rspace, ovl_bp_nspace, ovl_pct_nspace]
    if orient_match and ovl_pct_rspace > 50 and ovl_pct_nspace > 50:
        # not sure this could happen - these are quasi-identical alignments...
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'multimapper_self_id', overlaps)
        stats_record = ('issue', 'multimapper_self_id')
    elif not orient_match and ovl_pct_rspace > 50 and ovl_pct_nspace > 50:
        # again, quasi-identical alignments, but in opposite directions. For CN/ploidy > 1,
        # this could hold information...
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'multimapper_self_rev', overlaps)
        stats_record = ('issue', 'multimapper_self_rev')
        # this gets back into the Q since it has different orientation
        put_back = True
    elif ovl_pct_rspace > 50 and ovl_pct_nspace <= 50:
        # the same part of the read aligns to several non- or slightly overlapping
        # parts of a segment/node, potentially indicating a false duplication in the segment/node
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'false_duplication', overlaps)
        stats_record = ('issue', 'false_duplication')
    elif ovl_bp_rspace > 0 and ovl_pct_rspace < 50 and ovl_bp_nspace < 0:
        # a small-ish gap in the alignment could indicate a false expansion in the assembled
        # contig sequence
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'false_expansion', overlaps)
        stats_record = ('issue', 'false_expansion')
    elif ovl_pct_rspace <= 50 and ovl_pct_nspace > 50:
        # different parts of the read align all to the same part in the segment/node
        # that could be a sequence copy number error
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'false_copynumber', overlaps)
        stats_record = ('issue', 'false_copynumber')
    elif ovl_bp_rspace < 0 and ovl_bp_nspace > 0 and ovl_pct_nspace <= 50:
        # different parts of the read align to slightly overlapping parts of the node,
        # suggesting a sequence collapse
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'false_collapse', overlaps)
        stats_record = ('issue', 'false_collapse')
    elif ovl_bp_rspace > 0 and ovl_bp_nspace > 0:
        # this should be one single alignment, may still be caused by
        # assembly issues, but alignment or read error seems more plausible
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'alignment_error', overlaps)
        stats_record = ('issue', 'alignment_error')
    elif ovl_bp_rspace <= 0 and ovl_bp_nspace <= 0:
        # this should be one single alignment, but assembly error
        # in between the aligned blocks seems more likely. This is
        # in analogy to "broken_alignment_other"
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'broken_alignment_self', overlaps)
        stats_record = ('issue', 'broken_alignment_self')
    else:
        issue_record = IssueRecord.make_issue_record(record_last, record_first, 'unhandled', overlaps)
        raise RuntimeError(f'Unhandled self-alignment combination: {issue_record}')
    return issue_record, stats_record, put_back


def process_read_alignments(node_lengths, edge_lengths, orient_split, args):
    """
    Umbrella function to process all alignments of one single read
    (to the assembly graph) in one go.
    - Alignments are processed in order from left to right in read sequence coordinates
    --- for overlapping alignments, secondary sort by MAPQ in descending order
    - alignments to paths |P| > 1 count towards the support of the edges between
    each pair of connected nodes in the path
    --- sequence support is computed in that case, but is only approx. correct
    - if additional edges are induced by the long-read alignments, then these
    are only introduced between the current start node of the path and all last
    nodes of the last read alignment block
    --- special cases apply if there is an alignment disconnect
    """
    read_name, alignments = args

    sort_columns = ['read_align_start', 'read_align_end', 'mapq']
    sort_order = [True, True, False]

    # store read/node info
    # for last alignment block
    last_end_alignments = col.deque()

    # return list of node support
    node_support = []
    # return list of edge support
    edge_support = []
    # collect some stats on issues
    issue_records = []
    stats = col.Counter()

    get_nodes = fnt.partial(extract_path_members, orient_split)
    get_support = fnt.partial(collect_node_sequence_support, read_name, node_lengths, edge_lengths)

    for _, aln_record in alignments.sort_values(sort_columns, ascending=sort_order, inplace=False).iterrows():
        path_members = get_nodes(aln_record['path'])
        read_align_start = aln_record['read_align_start']
        read_align_end = aln_record['read_align_end']
        members_node_support, path_edge_support = get_support(
            path_members, aln_record
        )
        node_support.extend(members_node_support)
        edge_support.extend(path_edge_support)

        first_member = members_node_support[0]

        while 1:
            try:
                last_node, last_start, last_end = last_end_alignments.popleft()
            except IndexError:
                break
            if last_node is None:
                break

            # compute overlaps (ovl) between alignments
            # in "read space" (rspace) and in "node space" (nspace)
            # NB: alignments grouped by reads, i.e. it is always the same read
            # but not necessarily the same node/segment
            ovl_bp_rspace = min(read_align_end, last_end) - max(read_align_start, last_start)
            # if the overlap is negative (= alignments are apart), the percent overlap is 0
            ovl_pct_rspace = max(
                0,
                max(
                    round(ovl_bp_rspace / (read_align_end - read_align_start) * 100, 2),
                    round(ovl_bp_rspace / (last_end - last_start) * 100, 2)
                )
            )

            # need to record info if nodes are different
            same_node, same_orientation, ovl_bp_nspace, ovl_pct_nspace = check_intra_node_overlap(last_node, first_member)
            overlaps = [ovl_bp_rspace, ovl_pct_rspace, ovl_bp_nspace, ovl_pct_nspace]

            # a ton of case/switch...
            if ovl_bp_rspace < 0 and not same_node:
                # spaced alignments in read space between different nodes
                # >>> gap in assembly, induces an edge
                edge_specs = {
                    'length': ovl_bp_rspace,
                    'error_type': 'assembly_gap'
                }
                edge_support.append(Edge.make_edge(last_node, first_member, edge_specs))
                stats[('edge', 'assembly_gap')] += 1
            elif ovl_bp_rspace >= 0 and ovl_pct_rspace >= 50 and not same_node:
                # substantially overlapping alignments in read space between different nodes
                # >>> that is likely just a multimapper
                record_last = last_node.as_prefixed_record('a', last_start, last_end)
                record_first = first_member.as_prefixed_record('b', read_align_start, read_align_end)
                issue_records.append(IssueRecord.make_issue_record(record_last, record_first, 'multimapper_other', overlaps))
                stats[('issue', 'multimapper_other')] += 1
                # since this is a multimapper, the alignment needs to be considered
                # also for subsequent cases...
                last_end_alignments.append((last_node, last_start, last_end))
            elif ovl_bp_rspace >= 0 and ovl_pct_rspace < 50 and not same_node:
                # slightly overlapping alignments in read space between different nodes
                # >>> alignment broken (could be due to assembly errors at node boundaries)
                # make an edge
                edge_specs = {
                    'length': ovl_bp_rspace,
                    'error_type': 'broken_alignment_other'
                }
                edge_support.append(Edge.make_edge(last_node, first_member, edge_specs))
                stats[('edge', 'broken_alignment_other')] += 1
            elif same_node:
                # for now, assume all self(-looping) alignments are likely resulting from an assembly error,
                # or result from alignment problems in highly repetitive regions. At the moment,
                # these situation should not induce an edge (maybe consider later for error correction)
                issue_record, stats_record, put_back = handle_self_alignments(
                    last_node, last_start, last_end,
                    first_member, read_align_start, read_align_end,
                    same_orientation,
                    ovl_bp_rspace, ovl_pct_rspace,
                    ovl_bp_nspace, ovl_pct_nspace
                )
                if put_back:
                    last_end_alignments.append((last_node, last_start, last_end))
                issue_records.append(issue_record)
                stats[stats_record] += 1
            else:
                record_last = last_node.as_prefixed_record('a', last_start, last_end)
                record_first = first_member.as_prefixed_record('b', read_align_start, read_align_end)
                issue_record = IssueRecord.make_issue_record(record_last, record_first, 'unhandled', overlaps)
                raise RuntimeError(f'Unhandled alignment combination: {issue_record} ')
        # end of while block
        last_member = members_node_support[-1]
        last_end_alignments.append((last_member, read_align_start, read_align_end))
        last_end_alignments.append((None, -1, -1))

    return node_support, edge_support, issue_records, stats


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


def compute_gaf_record_id(row):
    row_key = ''.join(str(x) for x in row)
    row_hash = hashlib.md5(row_key.encode('ascii')).hexdigest()
    return row_hash


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
            raise ValueError(f'Illegal read alignment orientation: "{x}"')
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
    logger.debug('Adding row ID to GAF input...')
    alignments['_id'] = alignments.apply(compute_gaf_record_id, axis=1)
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
    logger.debug(f'Loaded {alignments.shape[0]} GAF records...')

    read_counter = 0
    for read_name, read_alignments in alignments.groupby('read_name'):
        read_counter += 1
        yield read_name, read_alignments
    logger.debug(f'Iterated through all alignments grouped by {read_counter} reads')
    return


def process_gaf_input(gaf_file_path, gaf_file_cache, ignore_cache, node_lengths, edge_lengths, num_jobs):
    """
    The current implementation only works if the GAF file has an informative MAPQ column
    """
    segment_re = re.compile(SEGMENT_ORIENTATION_SYMBOLS)
    process_read_aln = fnt.partial(process_read_alignments, node_lengths, edge_lengths, segment_re)

    complete_stats = col.Counter()
    edges = []
    nodes_support = []
    issues = []
    with mp.Pool(num_jobs) as pool:
        resit = pool.imap_unordered(process_read_aln, read_gaf_input(gaf_file_path, gaf_file_cache, ignore_cache))
        for node_support, edge_support, found_issues, stats in resit:
            complete_stats['reads'] += 1
            if complete_stats['reads'] % 10000 == 0:
                logger.info(f'Processed {complete_stats["reads"]} read alignment records')
            complete_stats['edges'] += len(edge_support)
            complete_stats['nodes_support'] += len(node_support)
            complete_stats['issues'] += len(found_issues)
            edges.extend(edge_support)
            nodes_support.extend(node_support)
            issues.extend(found_issues)
            complete_stats.update(stats)
    logger.debug('All GAF records processed')
    logger.debug(f'Total reads {complete_stats["reads"]} '
                 f'/ edges {complete_stats["edges"]} '
                 f'/ sequence support records {complete_stats["nodes_support"]}')

    return nodes_support, edges, issues, complete_stats


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
        '.node-support.cache.pck',
        '.edges-align.cache.pck',
    ]
    cache_names = [
        'gaf-input',
        'node-support',
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


def load_cached_graph_data(graph_cache):

    with open(graph_cache['node-lengths'], 'rb') as cache:
        node_lengths = pck.load(cache)
    logger.debug(f'Loaded cached length data for {len(node_lengths)} graph nodes (segments)')

    with open(graph_cache['edge-lengths'], 'rb') as cache:
        edge_lengths = pck.load(cache)
    logger.debug(f'Loaded cached edge length data for {len(edge_lengths)} graph edges (links)')

    return node_lengths, edge_lengths


def cache_alignment_data(node_support, edges, cache_files, dump_tables):

    logger.debug('Dumping alignment-derived edge cache...')
    with open(cache_files['edges-align'], 'wb') as cache:
        pck.dump(edges, cache)
    if dump_tables:
        logger.debug('Dumping alignment-derived edge cache as HDF/table...')
        hdf_path = cache_files['node-support'].with_suffix('.h5')
        with pd.HDFStore(hdf_path, 'w', complevel=9) as hdf:
            dump_df = pd.DataFrame.from_records(e.as_dict() for e in edges)
            hdf.put('cache', dump_df, format='fixed')

    logger.debug('Dumping node support cache...')
    with open(cache_files['node-support'], 'wb') as cache:
        pck.dump(node_support, cache)
    if dump_tables:
        logger.debug('Dumping node support cache as HDF/table...')
        hdf_path = cache_files['node-support'].with_suffix('.h5')
        with pd.HDFStore(hdf_path, 'w', complevel=9) as hdf:
            dump_df = pd.DataFrame.from_records(n.as_dict() for n in node_support)
            hdf.put('cache', dump_df, format='fixed')

    return


def dump_summary_statistics_file(stat_summary_path, cache_folder, cache_prefix, summary_stats):

    logger.debug('Dumping summary statistics file...')
    if stat_summary_path is None:
        stat_summary_path = cache_folder.parent / pl.Path(cache_prefix + '.stats.tsv')
        try:
            fd = open(stat_summary_path, 'w')
            fd.close()
            stat_summary_path.unlink()
        except PermissionError:
            logger.warning(
                f'The parent folder to the cache folder is not writable/accessible: {cache_folder.parent} '
                f'Creating the statistics summary output file in the cache folder: {cache_folder}'
            )
            stat_summary_path = cache_folder / pl.Path(cache_prefix + '.stats.tsv')
    stat_summary_path = stat_summary_path.resolve().absolute()
    logger.debug(f'Dumping summary statistics to path: {stat_summary_path}')
    stat_records = col.defaultdict(list)
    for k, v in summary_stats.items():
        try:
            item_type, item_info = k
        except ValueError:
            item_type = 'item'
            item_info = k
        stat_records[item_type].append((item_info, v))
    with open(stat_summary_path, 'w') as table:
        for item_type in ['item', 'edge', 'issue']:
            item_records = sorted(stat_records[item_type])
            for item_info, value in item_records:
                _ = table.write(f'{item_type}\t{item_info}\t{value}\n')
    logger.debug('Summary statistics saved.')
    return


def dump_issue_records_file(issue_records_path, cache_folder, cache_prefix, issue_records):

    logger.debug('Dumping issue records file...')
    if issue_records_path is None:
        issue_records_path = cache_folder.parent / pl.Path(cache_prefix + '.issues.tsv.gz')
        try:
            fd = open(issue_records_path, 'w')
            fd.close()
            issue_records_path.unlink()
        except PermissionError:
            logger.warning(
                f'The parent folder to the cache folder is not writable/accessible: {cache_folder.parent} '
                f'Creating the statistics summary output file in the cache folder: {cache_folder}'
            )
            issue_records_path = cache_folder / pl.Path(cache_prefix + '.issues.tsv.gz')
    issue_records_path = issue_records_path.resolve().absolute()
    logger.debug(f'Dumping issue records to path: {issue_records_path}')
    issue_table = pd.DataFrame.from_records(ir.as_dict() for ir in issue_records)
    issue_table.to_csv(issue_records_path, sep='\t', header=True, index=False)
    logger.debug('Issue records saved.')
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

        node_support, edges, issue_records, align_stats = process_gaf_input(
            args.align,
            cache_files['gaf-input'],
            args.ignore_cache,
            node_lengths,
            edge_lengths,
            args.num_jobs
        )
        dump_summary_statistics_file(args.stat_summary, args.cache_folder, args.cache_prefix, align_stats)
        dump_issue_records_file(args.issue_records, args.cache_folder, args.cache_prefix, issue_records)
        cache_alignment_data(node_support, edges, cache_files, args.cache_tables)
        logger.debug('Processing complete.')
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
