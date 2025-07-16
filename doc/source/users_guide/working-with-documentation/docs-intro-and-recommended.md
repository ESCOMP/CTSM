.. _docs-intro-and-recommended:

# Working with the CTSM documentation
.. _editing-the-documentation:

## One-time setup
You will need to have some software installed on your computer in order to build and the documentation and view the results:
- :ref:`building-docs-prereqs-mac`
- :ref:`building-docs-prereqs-windows`

## Editing the documentation
First, you will need a clone of CTSM to get all the documentation files and infrastructure. (If you're on Windows, you will make this clone in your :ref:`Ubuntu VM <install-wsl>`.) Note that you will clone this to your own computer, not Derecho or any cluster or anything.

The CTSM documentation is built from files in the `doc/source/tech_note/` and `doc/source/users_guide/` directories. These files are written in a mixture of what are called "markup languages." You may already be familiar—for better or for worse—with the LaTeX markup language. Fortunately, our documentation is simpler than that. It was originally written entirely in [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html), and it still mostly is, as you can tell by the predominance of .rst files. However, it's also possible to write Markdown documents (.md), which is nice because it's a much simpler and more widespread format (although see :ref:`tips-for-working-with-markdown`). If you've formatted text on GitHub, for instance, you've used Markdown.

Editing the documentation is as simple as opening the source file for the page you want to edit, then changing text. Make sure to use either reStructuredText or Markdown syntax, depending on the file's extension (.rst or .md, respectively). Note that "opening the source file" isn't completely straightforward on Windows; see :ref:`editing-text-files-wsl`.

If you're confident in your changes, or you're _not_ confident in your ability to preview and test the documentation (see [Building the documentation (recommended method)](#building-the-documentation) below), all you need to do is commit your changes and submit a pull request to the [CTSM GitHub repo](https://github.com/ESCOMP/CTSM). Automated testing will check the updated documentation for any errors, and a CTSM software engineer will review your PR. If everything looks good, they will merge it into the codebase and update the website.

.. _building-the-documentation:

## Building the documentation
We strongly suggest building the documentation on your personal computer before submitting a pull request, so that you can preview what your changes will look like. The recommended way to do this is using the `doc-builder` tool in conjunction with a "containerized" version of some required software.

### Directories
You will need a place to build the documentation. It's fine if that doesn't exist; the build tool will make it for you. The only restriction is that, at least for the recommended method described here, **your build directory must be somewhere in your CTSM clone**. (We recommend starting the name of your build directory with `_build` because CTSM knows to ignore such directories when it comes to `git`.) The instructions here assume you want to do your build in `doc/_build/`.

### Building the docs
All you need to do to build the docs with our recommended method is
```shell
cd doc
./build_docs -b _build -c -d
```

This runs a complicated series of scripts and software that culminate in something called Sphinx converting the .rst and .md files into HTML webpages.

The `-c` means "do a clean build." If you leave it off, Sphinx will only rebuild files it thinks have changed since the last time you built the docs. It's not always right, which can lead to problems. If you get unexpected errors without `-c`, rerunning with `-c` is the first troubleshooting step.

The `-d` means "run using the container." If you're not using the container and are instead using the `ctsm_pylib` Conda environment **(not recommended)**, leave off `-d`, and make sure you've activated `ctsm_pylib` before running the above command. See the "Container software or Conda environment" sections for :ref:`Mac <container-or-conda-mac>` or :ref:`Windows <container-or-conda-windows>` for more information on these two methods.

(Do `./build_docs --help` for more information and options.)

## Viewing your built docs

Note that there is a menu in the lower left of the webpage that lets readers switch between different versions of the documentation. The links to versions in this menu will not work when using the build command given above. If you wish to preview this version switching functionality, see :ref:`building-docs-multiple-versions`.

The process for viewing your build in a web browser differs depending on what kind of computer you have.

### Mac

You can open your build of the documentation in your default browser with
```shell
open _build/html/index.html
```

### Windows (Ubuntu VM)

Assuming you installed the WSL Utilities in the :ref:`windows-docs-ubuntu-utilities` setup step, you can open your build of the documentation like so:
```shell
wslview _build/html/index.html
```
If you didn't, you can do
```shell
explorer.exe $(wslpath -w _build/html/index.html)
```
These both do the same thing, but the `wslview` method is simpler. Either way, at least the first time you do this, it will open a window asking which app you'd like to view the HTML file in. Choose a browser like Microsoft Edge or Chrome. At the bottom of the window, you can then choose whether you always want to open HTML files using the selected app or just this once.
