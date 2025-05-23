.. _building-docs-multiple-versions:

Building multiple versions of the documentation
===============================================

There is a menu in the lower left of the webpage that lets readers switch between different versions of the documentation. Populating this menu involves a few steps.

First, look at the ``version_list`` line in ``docs/conf.py``. Edit that list as needed so it contains the name of each version you want.

Next, you will need to build the documentation once for each version. To build a version called ``latest``, you would first check out the corresponding version of the docs, then do:

.. code:: shell

   cd doc
   ./build_docs -r $HOME/path/to/build-dir -d -v latest

This will build the documentation in ``$HOME/path/to/build-dir/versions/latest``. Open ``$HOME/path/to/build-dir/versions/latest/html/index.html`` to see the result.

You can also leave off the ``-v latest``, in which case the current Git branch name will be used as the version name. Note, though, that in this case you will need to manually create the ``$HOME/path/to/build-dir/versions/branch-name/`` directory if it doesn’t already exist.
