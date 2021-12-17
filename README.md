# ONTourage
Home of the ONTourage tool

## Scope
`ONTourage` is in prototype stage. `ONTourage` is supposed to develop into a tool that supports users in exploring
genome graphs (now: HiFi assembly graphs) with matching ONT alignments.

For single chromosome graphs, `ONTourage` can easily run on standard laptops.

## Mode of operation
`ONTourage` caches all input data structures, i.e., the graph and the graph alignments. At the current stage of the
development, any version bump likely means you have to recreate the data cache.

## Setup
Clone this github repo, set up the Conda environment under `environment/dev_ontourage.yml`, activate the environment
and run `pip install --editable .` in the repo root folder. If everything worked out, you can run `ontourage --version`
from your shell and it should print the version number.

### Caching the graph / GFA

`ontourage -j NUM_CPU process-graph -g PATH_TO_GFA_FILE -cf PATH_TO_PUT_THE_CACHED_DATA -cp PREFIX_FOR_CACHE_FILE_NAMES`

You can force recreating the cache by appending the command `--ignore-cache` to the above command.
If you want to analyze the cached data with other tools, you can add the option `--cache-tables` to dump
the cached data in form of HDF5/h5 files (processable in Python using the `pandas` package).

### Caching the alignments / GAF
(has no use at the moment)

### Extracting (path) sequences from the graph

`ontourage extract-seq -cf PATH_TO_PUT_THE_CACHED_DATA -cp PREFIX_FOR_CACHE_FILE_NAMES -of FASTA_OUT_FILE -pn START,END -inv`

Traverse the (directed) graph after adding all reverse complement edges (parameter `-inv`), i.e., for an edge from A+ to B-,
also add the edge B+ to A- to the graph (happens on-the-fly). Adding the `-inv` option is highly recommended if the tool
producing the graph did not sanity-check the graph for (directed) edges potentially blocking traversals, e.g., by creating
nodes with only incoming edges in an otherwise linear path in a component of the graph.
The traversal starts at `START` and attempts to reach the node `END`. If the user does not specify any sequence
orientation (e.g., `START+` or `END-`), then all possible orientations for `START` and `END` are tested
(if they are part of the graph). If several paths/traversals are found, then the shortest one is dumped to FASTA.

In case no path can be found, or one wants to force the visit of certain nodes along the path, it is possible to specify
steps in the form of `START,STEP,STEP,...,END`. In that case, one can add the option `-skip` to allow skipping over steps
that cannot be resolved during traversal, and continue with the next step. For example, the path `A,B,C,D,E,F` may lack a
traversable path between `C,D`, and will --- using `-skip` --- then result in `A,B,C,E,F`.
