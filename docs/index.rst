Met-ID API
==========

Met-ID exposes a local HTTP API from its Rust/Axum server. The API provides
metabolite identification, database queries, imzML reading, and CCS prediction
through the configured prediction sidecar. The Python ``MetIDServer`` class is
a convenience client for starting the server and calling these routes.

This reference describes the behavior implemented by the supplied Rust and
Python code. Route names, accepted values, and response shapes below reflect
that code, including current limitations noted at the end.

Starting the server
-------------------

The Rust server binds to ``127.0.0.1`` on the port passed to ``serve_api``.
When started through the Python client, a free local port is selected unless
one is supplied explicitly. The client launches the release binary with
``--server --port <port>``, sets ``METID_CCS_SIDECAR`` to the prediction
sidecar executable, and waits for ``GET /warmup`` to return HTTP 200.

The default Python project directory is ``../src-tauri`` relative to the
client module. The client finds the Cargo package name in ``Cargo.toml`` to
locate ``target/release/<package-name>``. If the binary is missing or
``kill_conflicts=True``, it runs ``cargo build --release`` first.

Python quick start
------------------

The example below assumes the client module is importable as ``metid_api`` and
that its Python dependencies and native build tools are installed.

.. code-block:: python

   from metid_api import MetIDServer

   server = MetIDServer(port=8765)
   try:
       print(server.list_tables())
       results = server.predict_ccs(
           "CCO",
           adduct="[M+H]+",
       )
       print(results)
   finally:
       server.stop()

If no port is specified, the client selects an available local port. The
server's base URL is available as ``server.base_url``.

HTTP routes
-----------

All routes are served on the loopback interface. JSON request bodies should
use ``Content-Type: application/json``. The client sends JSON bodies on a few
``GET`` requests; see the implementation notes for compatibility details.

.. list-table:: Routes
   :header-rows: 1
   :widths: 18 18 64

   * - Method
     - Path
     - Purpose
   * - ``GET``
     - ``/``
     - Returns the text ``Hello from Tauri API!``.
   * - ``GET``
     - ``/warmup``
     - Health/readiness check. Returns ``OK``.
   * - ``GET``
     - ``/name``
     - Accepts a JSON string body and returns a greeting JSON string.
   * - ``POST``
     - ``/process``
     - Adds 5 to every number in a JSON array.
   * - ``POST``
     - ``/identify``
     - Identifies metabolites from m/z values and optional CCS values.
   * - ``GET``
     - ``/get_database``
     - Accepts a JSON string body naming ``tables`` or a table to fetch.
   * - ``POST``
     - ``/post_database``
     - Placeholder route; currently performs no database write.
   * - ``POST``
     - ``/cutom_sql``
     - Executes the supplied SQL statement. The route is spelled ``cutom`` in
       the current router.
   * - ``GET``
     - ``/search_by_name``
     - Accepts a JSON object body with ``name`` and ``ccs`` fields.
   * - ``POST``
     - ``/read_imzml``
     - Reads an imzML file and returns m/z and intensity arrays.
   * - ``POST``
     - ``/predict_ccs``
     - Predicts CCS and m/z for one or more SMILES/adduct pairs.

Health and greeting routes
--------------------------

``GET /``
~~~~~~~~~

Returns a plain-text greeting:

.. code-block:: text

   Hello from Tauri API!

``GET /warmup``
~~~~~~~~~~~~~~~

Returns HTTP 200 with the body ``OK`` when the server is ready. The Python
client polls this route during startup.

``GET /name``
~~~~~~~~~~~~~

The handler expects a JSON string in the request body and returns a JSON
string.

.. code-block:: http

   GET /name
   Content-Type: application/json

   "Patrik"

.. code-block:: json

   "Hello Patrik!"

``POST /process``
~~~~~~~~~~~~~~~~~

Accepts a JSON array of numbers and returns an array with 5 added to each
value.

.. code-block:: json

   [1.0, 2.5, -1.0]

.. code-block:: json

   [6.0, 7.5, 4.0]

Metabolite identification
-------------------------

``POST /identify`` accepts a JSON object. ``masses`` may be a number, numeric
string, or an array of numbers/numeric strings. ``mass_error`` and ``mzwindow``
are scalar numbers or numeric strings. ``typ`` and ``adducts`` accept a string
or an array of strings; numeric values are also accepted and converted to
strings by the parser.

.. code-block:: json

   {
     "met": "HMDB(All)",
     "mat": "FMP-10",
     "typ": "All",
     "adducts": "All",
     "mass_error": 0,
     "masses": [421.1911, 245.123],
     "ccss": [-1, -1],
     "mzwindow": 5
   }

Fields:

.. list-table:: ``/identify`` request fields
   :header-rows: 1
   :widths: 16 18 66

   * - Field
     - Required
     - Meaning and accepted values
   * - ``met``
     - Yes
     - Metabolome/database selector passed to the SQL handler (for example,
       ``HMDB(All)``).
   * - ``mat``
     - Yes
     - Matrix selector passed to the SQL handler (for example, ``FMP-10``).
   * - ``typ``
     - Yes
     - One type or a list of types. The special string ``All`` expands to
       ``Exogenous``, ``Endogenous``, and ``Unspecified``.
   * - ``adducts``
     - Yes
     - One value or a list of values. In this implementation the special
       string ``All`` maps to ``Phenolic Hydroxyls``.
   * - ``mass_error``
     - Yes
     - Scalar numeric value or numeric string, converted to a string before
       calling the SQL handler.
   * - ``masses``
     - Yes
     - One m/z value or a non-empty list of m/z values.
   * - ``ccss``
     - No
     - Optional CCS value(s). CCS is used only if the list length equals the
       m/z list length and its first value is not ``-1``. If omitted, empty,
       set to the ``-1`` sentinel, or length-mismatched, CCS is not used.
   * - ``mzwindow``
     - Yes
     - Scalar numeric value or numeric string passed as the m/z window.

The response is a JSON object whose properties are result-column names and
whose values are arrays of strings. Exact properties depend on the SQL
handler's returned identification fields. The Python ``identify`` helper
converts this object to an ``IdentificationDataFrame`` and sets ``names`` as
the DataFrame index.

Database routes
---------------

``GET /get_database``
~~~~~~~~~~~~~~~~~~~~

Accepts a JSON string body. Use ``"tables"`` to list table names, or pass a
table name to fetch its contents.

.. code-block:: http

   GET /get_database
   Content-Type: application/json

   "tables"

Table-list response:

.. code-block:: json

   {"Tables": ["metabolites", "matrices", "adducts"]}

For a table request, the response has this shape:

.. code-block:: json

   {
     "Table": {
       "col_names": ["name", "mz"],
       "rows": [["example", 123.4]]
     }
   }

SQLite values are represented as JSON nulls, numbers, strings, or base64
strings for BLOB values. The Python helpers ``list_metabolites``,
``list_matrices``, ``list_functional_groups``, and ``list_adducts`` convert
table results to pandas DataFrames.

``POST /post_database``
~~~~~~~~~~~~~~~~~~~~~~~

This route is currently a placeholder. It obtains a database connection but
does not add or modify records.

``POST /cutom_sql``
~~~~~~~~~~~~~~~~~~~

Accepts a JSON string containing SQL, executes it on the application
database, and returns an empty result wrapper. The route spelling is exactly
``/cutom_sql`` in the current server code. This is a powerful local operation;
send only trusted SQL and do not expose the server beyond loopback.

``GET /search_by_name``
~~~~~~~~~~~~~~~~~~~~~~

Accepts a JSON object body and searches metabolite names using a substring
match. Set ``ccs`` to ``true`` to include CCS columns in addition to m/z
columns.

.. code-block:: json

   {"name": "dopamine", "ccs": true}

The response uses the same ``Table`` shape as ``/get_database``. The Python
helper returns a pandas DataFrame indexed by ``name``.

imzML route
-----------

``POST /read_imzml`` accepts a file path and two numeric parameters:

.. code-block:: json

   {
     "path": "C:/data/sample.imzML",
     "n": 100,
     "tolerance": 0.01
   }

The response contains parallel arrays:

.. code-block:: json

   {
     "mz": [100.0, 101.0],
     "intensity": [1200.0, 450.0]
   }

The Python ``read_imzml`` helper returns these values as a DataFrame with
``mz`` and ``intensity`` columns. The path is read by the local server process,
so it must be accessible from that machine.

CCS prediction
--------------

``POST /predict_ccs`` accepts a single SMILES string or a list of SMILES
strings. The adduct property may be named ``adduct`` or ``adducts`` and may
be a string or a list of strings. If omitted, the default adduct is
``[M+H]+``.

.. code-block:: json

   {
     "smiles": ["CCO", "c1ccccc1"],
     "adducts": ["[M+H]+", "[M+Na]+"]
   }

One SMILES can be paired with several adducts, or several SMILES can share
one adduct. Otherwise, the two lists must have the same length. Empty SMILES
values are rejected. Blank adduct values use ``[M+H]+``.

.. code-block:: json

   [
     {
       "smiles": "CCO",
       "adduct": "[M+H]+",
       "predicted_ccs": 123.4,
       "mz": 47.049
     }
   ]

Each result includes the input SMILES, selected adduct, predicted CCS, and
calculated m/z. Predictions are sent to the sidecar sequentially.

The Python ``predict_ccs`` helper accepts a scalar or one-dimensional list of
SMILES and an optional scalar/list of adducts. It returns a pandas DataFrame.
The helper normalizes blank or unsupported adducts to ``[M+H]+`` before
sending the request. The Rust route itself reports invalid SMILES/adducts as
HTTP 400 and sidecar failures as HTTP 502.

Error responses
---------------

The identification API maps JSON and parsing errors to HTTP 400 and SQLite
errors to HTTP 500. The custom JSON error handler used by ``/process`` and
``/cutom_sql`` returns a JSON object with an ``error`` property: missing JSON
content type and invalid JSON data return HTTP 400; other JSON request
rejections return HTTP 422. CCS validation errors return HTTP 400, while
sidecar failures return HTTP 502.

The Python client raises ``MetIDApiError`` for HTTP failures when using its
shared error helper. It retries several operations before raising a
``RuntimeError``; CCS client errors in the 4xx range are converted to
``ValueError`` without retrying.

Python client reference
-----------------------

``MetIDServer(project_dir=None, port=None, binary_name=None,
kill_conflicts=True)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Starts a local API process. ``project_dir`` points to the Rust project,
``port`` selects the local listening port, ``binary_name`` overrides the
Cargo package name, and ``kill_conflicts`` controls whether the client tries
to stop a process already using the port and rebuilds the release binary.
The client starts a one-hour auto-shutdown timer by default.

Public helper methods:

* ``stop()`` — attempts ``POST /shutdown`` and then terminates the child
  process. The supplied Rust router does not currently register ``/shutdown``.
* ``identify(...)`` — calls ``/identify`` and returns an
  ``IdentificationDataFrame``.
* ``list_tables()`` — lists SQLite tables.
* ``list_metabolites()``, ``list_matrices()``,
  ``list_functional_groups()``, ``list_adducts()`` — return table data as
  pandas DataFrames.
* ``search_by_name(name, ccs=False)`` — searches names and returns a DataFrame.
* ``read_imzml(path, n, tolerance)`` — reads an imzML file into a DataFrame.
* ``predict_ccs(smiles, adduct=None)`` — predicts CCS and m/z into a DataFrame.
* ``return_name(name)`` and ``process_items(items)`` — call the simple example
  routes ``/name`` and ``/process``.

``IdentificationDataFrame.as_dict(mass_column="oMass")`` groups rows by the
requested mass column and returns a mapping from each mass to a dictionary
of identifications. Repeated or blank metabolite names receive generated
unique keys.

Implementation notes
--------------------

* ``/name``, ``/get_database``, and ``/search_by_name`` are registered as
  ``GET`` routes but expect JSON request bodies. Some HTTP clients, proxies,
  and servers do not handle GET bodies consistently. A future API revision
  should consider changing these to ``POST`` or query parameters.
* The client calls ``POST /shutdown``, but the supplied router has no such
  route; ``stop()`` still terminates its child process directly.
* The client code uses ``np.asarray`` in ``_as_1d_list`` but does not import
  NumPy in the supplied file. Also, ``_as_1d_list`` is defined inside the
  class but called as an unqualified function by ``identify``. Those issues
  must be corrected for the Python ``identify`` helper to work as written.
* The CCS retry path interpolates ``e`` after the exception handler. Python
  clears exception variables after an ``except`` block, so the final error
  path may raise ``NameError`` instead of the intended ``RuntimeError``.
* ``/search_by_name`` constructs SQL using the supplied name in a formatted
  query. Bind the value as a SQL parameter before exposing this route to
  untrusted input.
* ``/read_imzml`` unwraps numeric conversions and indexes returned pairs
  directly. Invalid numeric input or unexpected pair data may panic rather
  than return a structured client error.
* The ``_static`` directory must exist under ``docs/`` if
  ``html_static_path = ["_static"]`` remains in ``docs/conf.py``. Otherwise,
  Sphinx emits a warning.

