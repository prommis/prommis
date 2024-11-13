# Display flowsheet connectivity

This page describes a simple Python tool to calculate and display flowsheet (or model) connectivity.
The tool is available as a single module, `prommis.conn.connectivity`.
This module can be run from the [command-line interface](#CLI) or used programmatically through
its [Python API](#API).

The tool creates and process a connectivity matrix, where each cell at the intersection of an Arc and Unit.

The connectivity matrix format is:

|Arcs |Unit 1|Unit 2|...|Unit N|
|-----|------|------|---|------|
|Arc1 |-1    |0     |...|0     |
|Arc2 | 0    |1     |...|0     |
|...  | ...  |...   |...|...   |
|ArcN | 0    |1     |...|0     |

Where each cell at the intersection of an Arc (row i) and Unit (column j)
is either:
*  -1 meaning Arc(i) is an outlet of Unit(j),
*  1 meaning Arc(i) is an inlet for Unit(j),
*  0 meaning there is no connection

The tool can also produce output for the [Mermaid](https://mermaid.js.org/) diagramming tool,
which can produce visual representations of the flowsheet.

(CLI)=
## Command-line Interface

The command-line interface (CLI) is installed as a script named `connectivity`.
This script generates connectivity information from flowsheets,
or graphs of the flowsheet structure from connectivity information,
or both together -- i.e. graphs of the flowsheet structure from the flowsheet.

It has two main 'modes': csv and module.

The *csv* mode starts from a connectivity matrix stored in a file, and 
generates MermaidJS code (optionally wrapped by Markdown or HTML) that can be
rendered as a simple diagram of the flowsheet.

The *module* mode starts from a Python module, calls its 'build()' function to
build a flowsheet, then either (a) writes out the connectivity CSV for that flowsheet, or
(b) generates MermaidJS code as in the csv mode above. The second function
is a convenience and is equivalent to generating the CSV file then running again
in 'csv' mode with that file as input, which will be shown below.

You can explicitly indicate the mode with the --type/-t argument, though the
program will try to infer it as well (anything ending in ".csv" will be assumed to
be a CSV file, for example).

Example command-lines (showing the two modes):

```shell

# Generate the connectivity matrix in uky_conn.csv
connectivity prommis.uky.uky_flowsheet -O uky_conn.csv --to csv

# Generate the MermaidJS code wrapped in a HTML page that can be viewed in a
# browser without any further installation (MermaidJS is fetched from the network)
# The page will be called 'uky_conn.html' (since no filename was specified).
connectivity uky_conn.csv --to html

# Print the 'raw' MermaidJS code to the console instead of to a file
connectivity  uky_conn.csv --to mermaid --output-file "-"

# Print mermaid info to default file, with streams labeled
connectivity uky_conn.csv --to mermaid --labels
# (console)> Output in: uky_conn.mmd

```

For more information about MermaidJS, see http://mermaid.js.org

(API)=
## Python API

```{py:module} prommis.conn.connectivity
---
synopsis: Connectivity matrix for a flowsheet
---
```

### Functions

The two entry points for the API are `create_from_matrix` and `create_from_model`.
They both have similar behavior and inputs, and both will return an instance of the
class {py:class}`Connectivity` if there is no output file.

```{eval-rst}
.. autofunction:: create_from_matrix
```

```{eval-rst}
.. autofunction:: create_from_model
```

### Classes

```{eval-rst}
.. autoclass:: Connectivity
    :show-inheritance:
    :members: units, streams, connections
```

```{eval-rst}
.. autoclass:: ModelConnectivity
    :members:
```

```{eval-rst}
.. autoclass:: ConnectivityBuilder
    :members:
```

```{eval-rst}
.. autoclass:: Mermaid
    :members:
```
