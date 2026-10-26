Met-ID Python and R clients
===========================

The Met-ID clients let you use metabolite identification, CCS prediction, and
Met-ID database resources from Python and R. Each client starts or connects to
a local Met-ID server and returns results in data structures suited to its
language: pandas DataFrames in Python and data frames in R.

This guide focuses on programming workflows. You usually do not need to call
the HTTP routes yourself.

Before you start
----------------

You need the Met-ID desktop application and its matching CCS prediction
sidecar installed on the same computer where your Python or R session runs.
The clients start the local server process and wait until it is ready. The
server only needs to be reachable from the local machine.

Download the Python client (v>= 2.0.0) wheel from the `Met-ID GitHub releases page
<https://github.com/pbjarterot/Met-ID/releases>`_ and install it with pip:

.. code-block:: console

   python -m pip install path/to/metid_api.whl
   
The R client uses ``httr``, ``jsonlite``,
and ``processx``. The R client functions shown below must be available in your
session, either by loading the client package or sourcing its client code.

Python
------

Start and stop the local server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Import ``MetIDServer`` and use it as a context manager. The context manager
stops the child server process when the block ends, including when an error
occurs. By default, the client finds an installed Met-ID executable and the
matching ``prediction_model`` sidecar. It chooses a free local port and shuts
the server down after one hour of inactivity by default.

.. code-block:: python

   from metid_api import MetIDServer

   with MetIDServer() as metid:
       print(metid.base_url)
       print(metid.list_tables())

For explicit installation paths, pass ``server_path`` and ``sidecar_path`` or
set ``METID_SERVER_PATH`` and ``METID_CCS_SIDECAR`` in the environment. You can
also pass ``project_dir`` when the application is installed in a non-standard
location. ``port`` selects a specific local port; if that port is occupied,
the client raises an error. ``startup_timeout`` controls how long the client
waits for readiness, and ``auto_shutdown_after`` sets the idle shutdown timer
in seconds. Set ``auto_shutdown_after=None`` to disable that timer.

If you manage the process explicitly, call ``metid.stop()`` when finished:

.. code-block:: python

   metid = MetIDServer(port=8765)
   try:
       results = metid.identify(mz=[180.0634, 196.0583])
   finally:
       metid.stop()

Metabolite identification
~~~~~~~~~~~~~~~~~~~~~~~~

``identify`` accepts one or more m/z values and optional corresponding CCS
values. It returns an ``IdentificationDataFrame`` (a pandas DataFrame subclass)
with metabolite names as the index. The defaults select ``HMDB(All)`` and the
``FMP-10`` matrix.

.. code-block:: python

   with MetIDServer() as metid:
       hits = metid.identify(
           metabolome="HMDB(All)",
           matrix="FMP-10",
           typ="All",
           adducts="All",
           mass_error=5,
           mz=[180.0634, 196.0583],
           mzwindow=5,
       )

       print(hits.head())
       print(hits.columns)

Use ``ccs=-1`` (the default) or ``ccs=None`` to search without CCS. For a
multiple-mass query, provide one CCS value per m/z, or the single ``-1``
sentinel to disable CCS for every mass. Other length mismatches raise a
``ValueError``. The server only uses CCS filtering when the values are valid
and aligned with the m/z values.



The special value ``typ="All"`` asks the server to include Exogenous,
Endogenous, and Unspecified records. The current server maps ``adducts="All"``
to its configured default filter. You can also pass a string or list for these
selectors, depending on the database query you need.

To group identification rows by observed mass, use ``as_dict``. The default
grouping column is ``oMass``; pass another column name if your result uses a
different observed-mass field.


Identify output columns
------------------------
.. list-table:: Identification output columns
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``names``
     - Name of the annotated feature.
   * - ``oMass``
     - Observed m/z.
   * - ``aMass``
     - Adjusted m/z.
   * - ``tMass``
     - Theoretical m/z from the database.
   * - ``dMass``
     - Difference between adjusted and theoretical m/z.
   * - ``dPPM``
     - Difference between adjusted and theoretical m/z, in ppm.
   * - ``matrix``
     - Matrix used for identification.
   * - ``accession``
     - Database accession ID.
   * - ``smiles``
     - SMILES representation.
   * - ``msms``
     - Whether an MS/MS spectrum is available.
   * - ``formula``
     - Chemical formula.
   * - ``possible_derivs``
     - Number of annotated possible derivatizations.
   * - ``ccs``
     - Observed CCS.
   * - ``tCCS``
     - Theoretical CCS from the database.
   * - ``ccsPPM``
     - Difference between observed and theoretical CCS, in ppm.
   * - ``Posterior``
     - Posterior probability for the annotation.

.. code-block:: python

   grouped = hits.as_dict(mass_column="oMass")
   for mass, identifications in grouped.items():
       print(mass, list(identifications))

CCS prediction
~~~~~~~~~~~~~~

``predict_ccs`` accepts one SMILES string or a one-dimensional sequence. An
adduct can be omitted, supplied once for all SMILES, or supplied as a sequence
that pairs with the SMILES. A single SMILES can be paired with several adducts;
otherwise the two input lengths must match. The method returns a DataFrame
with ``smiles``, ``adduct``, ``predicted_ccs``, and ``mz`` columns.

.. code-block:: python

   with MetIDServer() as metid:
       predictions = metid.predict_ccs(
           ["CCO", "c1ccccc1"],
           adduct="[M+H]+",
       )
       print(predictions)

The Python wrapper accepts ``[M+H]+``, ``[M]+``, ``[M+Na]+``, and ``[M+K]+``.
Blank or unsupported adducts are replaced with the ``[M+H]+`` fallback before
the request is sent. Invalid SMILES are reported as a ``ValueError``.

Browse database data and search names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The client exposes common tables as DataFrames:

.. code-block:: python

   with MetIDServer() as metid:
       tables = metid.list_tables()
       metabolites = metid.list_metabolites()
       matrices = metid.list_matrices()
       groups = metid.list_functional_groups()
       adducts = metid.list_adducts()

       dopamine = metid.search_by_name("dopamine")
       print(dopamine.head())

``search_by_name(name, ccs=False)`` searches by name. Set ``ccs=True`` to ask
for CCS-related fields in addition to m/z fields. To retrieve another table,
use ``metid.list_tables()`` to find its name; the supplied Python client has
convenience methods for the tables listed above.

Read imzML data
~~~~~~~~~~~~~~~

``read_imzml(path, n, tolerance)`` reads an imzML file accessible to the local
Met-ID process and returns a DataFrame with ``mz`` and ``intensity`` columns.
The file path is interpreted on the computer running the server.

.. code-block:: python

   with MetIDServer() as metid:
       spectrum = metid.read_imzml(
           r"C:\data\sample.imzML",
           n=100,
           tolerance=0.01,
       )
       print(spectrum.head())

R
-

Configure and start the server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The R client starts the executable with ``processx``. Set ``METID_PROJECT_ROOT``
to the directory used by the R client's default path lookup. The default
lookup expects the release executable at
``met-id/src-tauri/target/release/Met-ID`` (``Met-ID.exe`` on Windows) and one
matching sidecar under ``met-id/src-tauri/pyinstaller/dist``.

For example, add the project location to ``~/.Renviron``:

.. code-block:: text

   METID_PROJECT_ROOT=/path/to/your/project

Alternatively, pass explicit ``binary`` and ``sidecar`` paths to
``start_metid_server``. The environment variables ``METID_SERVER_PATH`` and
``METID_CCS_SIDECAR_PATH`` provide the same overrides. The R client defaults
to port 8787; choose another port if it is already in use.

.. code-block:: r

   metid <- start_metid_server(port = 8787, timeout = 60)
   on.exit(stop_metid_server(metid), add = TRUE)

When running in a function, register ``on.exit`` immediately after startup so
the server is closed if later code fails. At the console, call
``stop_metid_server(metid)`` when you have finished.

Identify metabolites
~~~~~~~~~~~~~~~~~~~

``metid_identify`` returns a data frame. Use numeric vectors for multiple m/z
and CCS values. The default ``ccs = -1`` disables CCS filtering; ``ccs = NULL``
also disables it. For multiple m/z values, provide a matching CCS vector or
the scalar sentinel ``-1``.

.. code-block:: r

   metid <- start_metid_server()
   on.exit(stop_metid_server(metid), add = TRUE)

   hits <- metid_identify(
     metid,
     metabolome = "HMDB(All)",
     matrix = "FMP-10",
     typ = "All",
     adducts = "All",
     mass_error = 5,
     mz = c(180.0634, 196.0583),
     mzwindow = 5,
     ccs = c(-1, -1)
   )

   head(hits)

The server's ``typ="All"`` setting includes Exogenous, Endogenous, and
Unspecified records. In the current server implementation,
``adducts="All"`` maps to its configured default filter.

Predict CCS
~~~~~~~~~~~

``metid_predict_ccs`` returns one row per prediction, with the input SMILES,
adduct, predicted CCS, and calculated m/z. If ``adduct`` is omitted, the
server uses ``[M+H]+``.

.. code-block:: r

   predictions <- metid_predict_ccs(
     metid,
     smiles = c("CCO", "c1ccccc1"),
     adduct = c("[M+H]+", "[M+H]+")
   )

   predictions

A single adduct can be used for multiple SMILES. For multiple adducts per
molecule, repeat the SMILES value for each requested pairing. Use adduct values
supported by the installed sidecar.

Work with database tables and name search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``metid_list_tables`` to inspect the database tables available in the
installed application. Common convenience functions return data frames:

.. code-block:: r

   tables <- metid_list_tables(metid)
   metabolites <- metid_list_metabolites(metid)
   matrices <- metid_list_matrices(metid)
   groups <- metid_list_functional_groups(metid)
   adducts <- metid_list_adducts(metid)

   dopamine <- metid_search_by_name(metid, "dopamine", ccs = TRUE)

For a table without a convenience function, call ``metid_get_table(metid,
"table_name")``. ``metid_search_by_name`` accepts a name and an optional
``ccs`` flag.

Read imzML data
~~~~~~~~~~~~~~~

``metid_read_imzml`` returns a data frame with parallel ``mz`` and
``intensity`` columns. The path must be visible to the local server process.

.. code-block:: r

   spectrum <- metid_read_imzml(
     metid,
     path = "data/sample.imzML",
     n = 100,
     tolerance = 0.01
   )
   head(spectrum)

Common client operations
------------------------

| Task | Python | R |
| --- | --- | --- |
| Start server | ``MetIDServer()`` | ``start_metid_server()`` |
| Stop server | ``stop()`` or context manager | ``stop_metid_server(client)`` |
| Identify | ``identify(...)`` | ``metid_identify(client, ...)`` |
| Predict CCS | ``predict_ccs(...)`` | ``metid_predict_ccs(client, ...)`` |
| List tables | ``list_tables()`` | ``metid_list_tables(client)`` |
| Read a table | convenience methods | ``metid_get_table(client, name)`` |
| Search names | ``search_by_name(...)`` | ``metid_search_by_name(client, ...)`` |
| Read imzML | ``read_imzml(...)`` | ``metid_read_imzml(client, ...)`` |

Errors and troubleshooting
--------------------------

* **Executable not found:** pass ``server_path`` and ``sidecar_path`` in
  Python, or ``binary`` and ``sidecar`` in R. Check the corresponding
  environment variables and executable permissions.
* **Server fails to start:** verify the executable and sidecar came from the
  same Met-ID release. In R, confirm ``METID_PROJECT_ROOT`` points to the
  expected directory layout. In Python, inspect the startup output printed by
  the client.
* **Port is already in use:** choose another ``port``. The clients do not need
  a public or network-facing port.
* **CCS prediction fails:** check that each SMILES is valid and the adduct is
  supported by the installed sidecar. Python normalizes unsupported adducts to
  ``[M+H]+``; R passes the selected adduct through to the server.
* **imzML file cannot be read:** use an existing path readable by the local
  server process, not merely by a different machine or remote R/Python client.
* **Identification returns no rows:** verify the selected metabolome, matrix,
  mass error, m/z window, and whether CCS filtering is enabled. Use ``-1`` to
  disable CCS filtering.

The clients retry some transient requests and report server errors as
exceptions. Python exposes ``MetIDApiError`` for HTTP failures; the CCS helper
raises ``ValueError`` for invalid client input. R reports request failures as
``metid_api_error`` conditions, which include a ``status_code`` when an HTTP
status was returned.

Low-level HTTP access
--------------------

The clients call the local Axum API internally. Most application code should
use the language-level functions above. Advanced users who need to implement
another client can inspect the server's route definitions: ``/identify`` for
metabolite search, ``/predict_ccs`` for CCS prediction, ``/get_database`` for
table access, ``/search_by_name`` for name lookup, and ``/read_imzml`` for
imzML input. These routes are local implementation details; the Python and R
client methods provide the recommended interface.
