.. _installPypi:

=========================
Installation via PyPI
=========================


pandas can be installed via pip from
`PyPI <https://pypi.org/project/pandas>`__.

::

    pip install pandas

Installing with ActivePython
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation instructions for
`ActivePython <https://www.activestate.com/activepython>`__ can be found
`here <https://www.activestate.com/activepython/downloads>`__. Versions
2.7, 3.5 and 3.6 include pandas.

Installing using your Linux distribution's package manager.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The commands in this table will install pandas for Python 3 from your distribution.

.. csv-table::
    :header: "Distribution", "Status", "Download / Repository Link", "Install method"
    :widths: 10, 10, 20, 50


    Debian, stable, `official Debian repository <https://packages.debian.org/search?keywords=pandas&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install python3-pandas``
    Debian & Ubuntu, unstable (latest packages), `NeuroDebian <http://neuro.debian.net/index.html#how-to-use-this-repository>`__ , ``sudo apt-get install python3-pandas``
    Ubuntu, stable, `official Ubuntu repository <https://packages.ubuntu.com/search?keywords=pandas&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install python3-pandas``
    OpenSuse, stable, `OpenSuse Repository  <https://software.opensuse.org/package/python-pandas?search_term=pandas>`__ , ``zypper in python3-pandas``
    Fedora, stable, `official Fedora repository  <https://admin.fedoraproject.org/pkgdb/package/rpms/python-pandas/>`__ , ``dnf install python3-pandas``
    Centos/RHEL, stable, `EPEL repository <https://admin.fedoraproject.org/pkgdb/package/rpms/python-pandas/>`__ , ``yum install python3-pandas``

**However**, the packages in the linux package managers are often a few versions behind, so
to get the newest version of pandas, it's recommended to install using the ``pip`` or ``conda``
methods described above.

Handling ImportErrors
~~~~~~~~~~~~~~~~~~~~~~

If you encounter an ImportError, it usually means that Python couldn't find pandas in the list of available
libraries. Python internally has a list of directories it searches through, to find packages. You can
obtain these directories with::

            import sys
            sys.path

One way you could be encountering this error is if you have multiple Python installations on your system
and you don't have pandas installed in the Python installation you're currently using.
In Linux/Mac you can run ``which python`` on your terminal and it will tell you which Python installation you're
using. If it's something like "/usr/bin/python", you're using the Python from the system, which is not recommended.

It is highly recommended to use ``conda``, for quick installation and for package and dependency updates.
You can find simple installation instructions for pandas in this document: ``installation instructions </getting_started.html>``.
