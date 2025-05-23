.. _docs-intro-and-recommended:

# Working with the CTSM documentation
.. _editing-the-documentation:

## One-time setup
In addition to a CTSM checkout in which to work, you will need to have some software installed in order to build the documentation:
- :ref:`building-docs-prereqs-mac`
- :ref:`building-docs-prereqs-windows`

## Editing the documentation
The CTSM documentation is built from files in the `doc/source/tech_note/` and `doc/source/users_guide/` directories. These files are written in a mixture of what are called "markup languages." You may already be familiar—for better or for worse—with the LaTeX markup language. Fortunately, our documentation is simpler than that. It was originally written entirely in [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html), and it still mostly is, as you can tell by the predominance of .rst files. However, it's also possible to write Markdown documents (.md), which is nice because it's a much simpler and more widespread format (although see :ref:`tips-for-working-with-markdown`). If you've formatted text on GitHub, for instance, you've used Markdown.

Editing the documentation is as simple as opening the source file for the page you want to edit, then changing text. Make sure to use either reStructuredText or Markdown syntax, depending on the file's extension (.rst or .md, respectively).

If you're confident in your changes, or you're _not_ confident in your ability to preview and test the documentation (see [Building the documentation (recommended method)](#building-the-documentation-recommended-method) below), all you need to do is commit your changes and submit a pull request to the [CTSM GitHub repo](https://github.com/ESCOMP/CTSM). Automated testing will check the updated documentation for any errors, and a CTSM software engineer will review your PR. If everything looks good, they will merge it into the codebase and update the website.

.. _building-the-documentation-recommended-method:

## Building the documentation (recommended method)
We strongly suggest building the documentation on your personal computer before submitting a pull request, so that you can preview what your changes will look like. The recommended way to do this is using the `doc-builder` tool in conjunction with a "containerized" version of some required software.

### Directories
You will need a place to build the documentation. It's fine if that doesn't exist; the build tool will make it for you. The only restriction is that, at least for the recommended method described here, **your build directory must be somewhere in your user home directory**, which we represent as `$HOME`. The instructions here assume you want to do your build in `$HOME/path/to/build-dir/`.

Your CTSM clone, from which you're building the documentation, also needs to be somewhere in your user home directory.

### Building a preview
Ensure that Docker Desktop is running. (You do not actually need to do anything in Docker Desktop; it just needs to be running. It's fine for it to be minimized or hidden.) Then all you need to do is
```shell
cd doc
./build_docs -b $HOME/path/to/build-dir -c -d
```

(Do `./build_docs --help` for more information and options.)

You can then open the documentation in a web browser by browsing to `$HOME/path/to/build-dir/html/` and opening `index.html`.

Note that there is a menu in the lower left of the webpage that lets readers switch between different versions of the documentation. The links to versions in this menu will not work when using the build command given above. If you wish to preview this version switching functionality, or you're building the docs in the process of actually updating the website, see :ref:`building-docs-multiple-versions`.
