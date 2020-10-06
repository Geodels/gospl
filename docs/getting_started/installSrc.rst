.. _installSrc:

=========================
Installation via Source
=========================

Python version support
----------------------

Officially Python 3.7.1 and above, 3.8, and 3.9.


See the :ref:`contributing guide <contributing>` for complete instructions on building from the git source tree. Further, see :ref:`creating a development environment <contributing.dev_env>` if you wish to create a *pandas* development environment.

Running the test suite
----------------------

pandas is equipped with an exhaustive set of unit tests, covering about 97% of
the code base as of this writing. To run it on your machine to verify that
everything is working (and that you have all of the dependencies, soft and hard,
installed), make sure you have `pytest
<https://docs.pytest.org/en/latest/>`__ >= 5.0.1 and `Hypothesis
<https://hypothesis.readthedocs.io/>`__ >= 3.58, then run:

::

    >>> pd.test()
    running: pytest --skip-slow --skip-network C:\Users\TP\Anaconda3\envs\py36\lib\site-packages\pandas
    ============================= test session starts =============================
    platform win32 -- Python 3.6.2, pytest-3.6.0, py-1.4.34, pluggy-0.4.0
    rootdir: C:\Users\TP\Documents\Python\pandasdev\pandas, inifile: setup.cfg
    collected 12145 items / 3 skipped

    ..................................................................S......
    ........S................................................................
    .........................................................................

    ==================== 12130 passed, 12 skipped in 368.339 seconds =====================

.. _install.dependencies:

Dependencies
------------

================================================================ ==========================
Package                                                          Minimum supported version
================================================================ ==========================
`setuptools <https://setuptools.readthedocs.io/en/latest/>`__    24.2.0
`NumPy <https://www.numpy.org>`__                                1.16.5
`python-dateutil <https://dateutil.readthedocs.io/en/stable/>`__ 2.7.3
`pytz <https://pypi.org/project/pytz/>`__                        2017.3
================================================================ ==========================

.. _install.recommended_dependencies:

Recommended dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

* `numexpr <https://github.com/pydata/numexpr>`__: for accelerating certain numerical operations.
  ``numexpr`` uses multiple cores as well as smart chunking and caching to achieve large speedups.
  If installed, must be Version 2.6.8 or higher.

* `bottleneck <https://github.com/pydata/bottleneck>`__: for accelerating certain types of ``nan``
  evaluations. ``bottleneck`` uses specialized cython routines to achieve large speedups. If installed,
  must be Version 1.2.1 or higher.

.. note::

   You are highly encouraged to install these libraries, as they provide speed improvements, especially
   when working with large data sets.


.. _install.optional_dependencies:

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

Pandas has many optional dependencies that are only used for specific methods.
For example, :func:`pandas.read_hdf` requires the ``pytables`` package, while
:meth:`DataFrame.to_markdown` requires the ``tabulate`` package. If the
optional dependency is not installed, pandas will raise an ``ImportError`` when
the method requiring that dependency is called.

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
BeautifulSoup4            4.6.0              HTML parser for read_html (see :ref:`note <optional_html>`)
Jinja2                    2.10               Conditional formatting with DataFrame.style
PyQt4                                        Clipboard I/O
PyQt5                                        Clipboard I/O
PyTables                  3.4.4              HDF5-based reading / writing
SQLAlchemy                1.2.8              SQL support for databases other than sqlite
SciPy                     1.12.0             Miscellaneous statistical functions
xlsxwriter                1.0.2              Excel writing
blosc                     1.14.3             Compression for HDF5
fsspec                    0.7.4              Handling files aside from local and HTTP
fastparquet               0.3.2              Parquet reading / writing
gcsfs                     0.6.0              Google Cloud Storage access
html5lib                  1.0.1              HTML parser for read_html (see :ref:`note <optional_html>`)
lxml                      4.3.0              HTML parser for read_html (see :ref:`note <optional_html>`)
matplotlib                2.2.3              Visualization
numba                     0.46.0             Alternative execution engine for rolling operations
openpyxl                  2.6.0              Reading / writing for xlsx files
pandas-gbq                0.12.0             Google Big Query access
psycopg2                  2.7                PostgreSQL engine for sqlalchemy
pyarrow                   0.15.0             Parquet, ORC, and feather reading / writing
pymysql                   0.7.11             MySQL engine for sqlalchemy
pyreadstat                                   SPSS files (.sav) reading
pytables                  3.4.4              HDF5 reading / writing
pyxlsb                    1.0.6              Reading for xlsb files
qtpy                                         Clipboard I/O
s3fs                      0.4.0              Amazon S3 access
tabulate                  0.8.3              Printing in Markdown-friendly format (see `tabulate`_)
xarray                    0.12.0             pandas-like API for N-dimensional data
xclip                                        Clipboard I/O on linux
xlrd                      1.2.0              Excel reading
xlwt                      1.3.0              Excel writing
xsel                                         Clipboard I/O on linux
zlib                                         Compression for HDF5
========================= ================== =============================================================

.. _optional_html:

Optional dependencies for parsing HTML
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the following combinations of libraries is needed to use the
top-level :func:`~pandas.read_html` function:

* `BeautifulSoup4`_ and `html5lib`_
* `BeautifulSoup4`_ and `lxml`_
* `BeautifulSoup4`_ and `html5lib`_ and `lxml`_
* Only `lxml`_, although see :ref:`HTML Table Parsing <io.html.gotchas>`
  for reasons as to why you should probably **not** take this approach.

.. warning::

    * if you install `BeautifulSoup4`_ you must install either
      `lxml`_ or `html5lib`_ or both.
      :func:`~pandas.read_html` will **not** work with *only*
      `BeautifulSoup4`_ installed.
    * You are highly encouraged to read :ref:`HTML Table Parsing gotchas <io.html.gotchas>`.
      It explains issues surrounding the installation and
      usage of the above three libraries.

.. _html5lib: https://github.com/html5lib/html5lib-python
.. _BeautifulSoup4: https://www.crummy.com/software/BeautifulSoup
.. _lxml: https://lxml.de
.. _tabulate: https://github.com/astanin/python-tabulate
