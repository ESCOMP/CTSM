.. _building-docs-multiple-versions:

Building multiple versions of the documentation
===============================================

There is a menu in the lower left of the webpage that lets readers switch between different versions of the documentation. To build a website with this menu properly set up—so that all our versions appear and all the links work—you need to use ``docs/build_docs_to_publish`` instead of ``docs/build_docs``.

Note that this is not necessary in order for you to contribute an update to the documentation. GitHub will test this automatically when you open a PR. But if you'd like to try, this will generate a local site for you in ``_publish/`` and then open it:

.. literalinclude:: ../../../test/test_container_eq_ctsm_pylib.sh
   :start-at: ./build_docs_to_publish
   :end-before: VERSION LINKS WILL NOT RESOLVE
   :append: CMD _publish/index.html  # where CMD is open for Mac or wslview for Windows (Ubuntu VM)

**Note:** This is not yet supported with Podman on Linux (including Ubuntu VM on Windows). See `doc-builder Issue #27: build_docs_to_publish fails on Linux (maybe just Ubuntu?) with Podman <https://github.com/ESMCI/doc-builder/issues/27>`_. It does work with Docker on Linux, though.


How this works
--------------

``build_docs_to_publish`` loops through the ``VERSION_LIST`` variable in ``doc/version_list.py``:

.. literalinclude:: ../../../version_list.py
   :start-at: version of certain files we want to preserve
   :end-before: End version definitions

For each member of ``VERSION_LIST``, ``build_docs_to_publish`` checks out its ``ref``, then builds the documentation in a build directory. (``LATEST_REF`` is set because some files, folders, and submodules are important for how the build works and need to stay the same for each build.) Once the build is complete, ``build_docs_to_publish`` should reset your local repo copy (CTSM clone) to how it was before you called ``build_docs_to_publish``.

Next, ``build_docs_to_publish`` moves the HTML files from the build directory to the publish directory. The publish directory has a structure that matches the paths in the version dropdown menu's links. If a member of ``VERSION_LIST`` has ``landing_version=True``, its HTML will be at the top level. That makes it simple for people to find the default version of the docs at https://escomp.github.io/CTSM, rather than having to drill down further into something like ``https://escomp.github.io/CTSM/versions/latest``.
