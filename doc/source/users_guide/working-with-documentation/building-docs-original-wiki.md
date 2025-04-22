.. _building-docs-original-wiki:

# ⚠️ Original docs documentation from the GitHub Wiki

.. todo::

   ⚠️⚠️⚠️WARNING⚠️⚠️⚠️
   This page contains documentation that (a) is more complicated than you probably require and (b) has not been fully checked for accuracy with the latest documentation setup. Unless you have a very good reason, you should probably go to :ref:`docs-intro-and-recommended`.


## Table of contents

* [Intro to the CTSM documentation](#intro-to-the-ctsm-documentation)
* [What do I, as a contributor, need to do to update the tech note for my new feature?](#what-do-i-as-a-contributor-need-to-do-to-update-the-tech-note-for-my-new-feature)
* [Quick start to building the documentation](#quick-start-to-building-the-documentation)
* [Overview of the recommended build method and alternative methods](#overview-of-the-recommended-build-method-and-alternative-methods)
* [One-time setup needed for a given machine](#one-time-setup-needed-for-a-given-machine)
   * [Prerequisites](#prerequisites)
   * [Install Docker and download the required container](#install-docker-and-download-the-required-container)
   * [Install git-lfs](#install-git-lfs)
* [Management of image files](#management-of-image-files)
   * [Obtaining the image files](#obtaining-the-image-files)
   * [Seeing what files are tracked by git-lfs](#seeing-what-files-are-tracked-by-git-lfs)
   * [Adding a new image file type](#adding-a-new-image-file-type)
* [Procedure for building the html documentation](#procedure-for-building-the-html-documentation)
   * [Recommended directory structure for these directions](#recommended-directory-structure-for-these-directions)
   * [Initial steps for building the documentation](#initial-steps-for-building-the-documentation)
   * [Building a test version of the documentation for your own review](#building-a-test-version-of-the-documentation-for-your-own-review)
   * [Previewing the built documentation](#previewing-the-built-documentation)
   * [Updating the official html documentation](#updating-the-official-html-documentation)
   * [Dealing with errors](#dealing-with-errors)
      * [Input/output error](#inputoutput-error)
   * [Adding a new version](#adding-a-new-version)
* [Building a pdf of the tech note](#building-a-pdf-of-the-tech-note)
* [Resources for learning markup with reStructuredText and using Sphinx.](#resources-for-learning-markup-with-restructuredtext-and-using-sphinx)
* [Appendix: Other build methods](#appendix-other-build-methods)
   * [Running build_docs from all of these methods](#running-build_docs-from-all-of-these-methods)
   * [Relying more heavily on the Docker container](#relying-more-heavily-on-the-docker-container)
      * [Launching the Docker image](#launching-the-docker-image)
      * [Building the documentation](#building-the-documentation)
      * [Viewing the built documentation](#viewing-the-built-documentation)
      * [Committing to git repositories](#committing-to-git-repositories)
   * [Installing required software on your desktop/laptop (instructions for Mac)](#installing-required-software-on-your-desktoplaptop-instructions-for-mac)
      * [Prerequisites / assumptions](#prerequisites--assumptions)
      * [Installing Sphinx and the necessary Sphinx theme](#installing-sphinx-and-the-necessary-sphinx-theme)
      * [Installing latexmk](#installing-latexmk)
      * [Optionally installing components needed to build the PDF](#optionally-installing-components-needed-to-build-the-pdf)
   * [Building the documentation on cheyenne](#building-the-documentation-on-cheyenne)
* [Appendix: Editing tips](#appendix-editing-tips)

.. _intro-to-the-ctsm-documentation:

## Intro to the CTSM documentation

The CTSM documentation is written using [reStructuredText markup language](http://www.sphinx-doc.org/en/stable/rest.html). ReStructuredText is a markup language, like HTML, markdown or latex. (See below for more resources on learning reStructuredText.)

[Sphinx](http://www.sphinx-doc.org/en/stable/index.html) is a python tool for publishing reStructuredText documents in other formats such as HTML and PDF.

The CTSM documentation source is stored in the `doc/source` directory of the main CTSM repository. The built documentation is stored in the gh-pages branch of the [ESCOMP/ctsm-docs repository](https://github.com/escomp/ctsm-docs). The Sphinx-generated HTML pages are accessible from URL [https://escomp.github.io/ctsm-docs/](https://escomp.github.io/ctsm-docs/)

.. _what-do-i-as-a-contributor-need-to-do-to-update-the-tech-note-for-my-new-feature:

## What do I, as a contributor, need to do to update the tech note for my new feature?

As a contributor to CTSM, you don't necessarily need to know all of the information laid out below. If you need to update the tech note for some changes you have made, you should add or edit the necessary rst files in `doc/source/tech_note` on your CTSM branch. Then simply commit your changes to these rst files just as you would do for changes to Fortran source files, and include these rst changes in your Pull Request. Ideally these documentation changes will be on the same branch and the same Pull Request as the related source file changes, but they can also be made later, on a separate branch, and submitted to us via a separate Pull Request. In either case, the CTSM maintainers will review the documentation changes in the Pull Request and merge these changes to CTSM's master branch when ready.

You do *not* need to open a Pull Request for the built documentation: one of the CTSM maintainers will rebuild the documentation later. However, if you have made more than minor edits to the documentation, we prefer if you have tested the documentation build to make sure there are no errors in the edited rst files and that the html documentation appears how you intended. But you do *not* need to send us or point us to your rebuilt documentation: One of the CTSM maintainers will build the documentation themselves to review the build.

If you need to build the documentation, and/or if you need to change or add any images, you will need to keep reading for instructions on installing and using `git lfs`.

.. _quick-start-to-building-the-documentation:

## Quick start to building the documentation

This documentation gives detailed explanations of a number of workflows that can be used in building the documentation. This section briefly summarizes a single, recommended workflow, so that you can get started quickly.

First, as one-time setup, you will need to install:
- A recent version of git
- A recent version of python (python3.5 or later), available as `python3`
- [Docker](https://www.docker.com/products/docker-desktop)
- [Git LFS](https://git-lfs.github.com/), including running `git lfs install` (required if you will be adding or editing any image files, otherwise you can skip Git LFS for now) You should then launch Docker's desktop application, and run `docker pull escomp/base`.

Then (if you haven't already done so), clone the CTSM repository (this example assumes that it appears in a path like this: `~/ctsm-repos/ctsm`). Check out the branch from which you want to build the documentation, then run `./bin/git-fleximod update --optional`.

Then build the documentation like this:

```shell
mkdir ~/ctsm-repos/ctsm-docs
cd ~/ctsm-repos/ctsm/doc
./build_docs -b ~/ctsm-repos/ctsm-docs -d
```

and view it by opening the file `~/ctsm-repos/ctsm-docs/html/index.html`.

.. _overview-of-the-recommended-build-method-and-alternative-methods:

## Overview of the recommended build method and alternative methods

The directions here assume use of a documentation build method where you are working on your local desktop / laptop, and have some basic software installed there, but use a Docker container-based method for doing the actual documentation build. Other methods are documented in an [Appendix below](#appendix-other-build-methods). We prefer the Docker-based method because it avoids the complexity of installing multiple packages on your local machine, and ensures that all people building the documentation are using the same versions of Sphinx and other utilities. The main method described here leverages this Docker container while still allowing you to remain in the comfort of your own local environment. However, some people may prefer [the more Docker-centric method described in the Appendix](#relying-more-heavily-on-the-docker-container), because it has even fewer initial installation requirements. **That solution may also work more smoothly on a Windows system, though we have run the primary method successfully on Windows after ensuring that we have an executable named python3 installed.**

.. _one-time-setup-needed-for-a-given-machine:

## One-time setup needed for a given machine

.. _prerequisites:

### Prerequisites

These instructions assume that you have the following available on your machine:

- A recent version of git, ideally configured according to the recommendations in <https://github.com/ESCOMP/CTSM/wiki/Recommended-git-setup#one-time-git-configuration>.

- A recent version of python (python3.5 or later)
  - Note that this must be available as `python3` (not just `python`). If your python installation did not create a `python3` command, you will need to create an alias or symbolic link so that python can be invoked via `python3`. On Windows, one developer reports success after copying `C:\Users\USERNAME\AppData\Local\Programs\Python\Python39\python.exe` to `python3.exe` (in the same location).

.. _install-docker-and-download-the-required-container:

### Install Docker and download the required container

Download Docker to your personal desktop or laptop from [Docker](https://www.docker.com/products/docker-desktop) and follow the instructions. Then launch the Docker application.

You can then obtain the necessary Docker container by running:

```shell
docker pull escomp/base
```

Note: For some versions of Docker on a Mac, Docker's CPU usage can remain at 100% even after the documentation build process exits. This seems to be connected with the issue reported in <https://github.com/docker/for-mac/issues/5164>. The issue probably arises because the `build_docs` script mounts your entire home directory in the docker container. You can work around this issue by quitting Docker when you are done building the docs, or by unchecking the option, "Use gRPC FUSE for file sharing" in Docker's preferences. However, note that unchecking that option might have some negative consequences, so you may want to re-enable it when you are done building the documentation if you use Docker for other purposes as well.

.. _install-git-lfs:

### Install git-lfs

Git-LFS (Large File Support) is needed if you will be adding or editing any image files. **If there is even a chance that you will be doing so, it is best to go ahead and install this now, while you are thinking about it. This will prevent you from accidentally committing an image file directly to the repository, which is something we try very hard to avoid.** This is also needed for building the documentation if *not* using the Docker-based method.

Follow the instructions on the [Git LFS page](https://git-lfs.github.com/) for installing git-lfs on your platform. (On a Mac using homebrew, this can be done with `brew install git-lfs`.)

**Although the actual installation only needs to be done once per machine, each user who wants to use it must then run the following command to set up git-lfs for your user account on this machine:**

```shell
git lfs install
```

.. _management-of-image-files:

## Management of image files

**If you are adding or changing any image files, be sure you have installed git-lfs according to the [above instructions](#Install-git-lfs), including running** `git lfs install`. If you are not sure whether you have already run `git lfs install`, it is safe to rerun that command to be sure.

Image files are tracked using Git LFS (Large File Support). The images are stored somewhere on GitHub, but aren't actually part of the main repository. Instead, the repository contains small text files that direct git-lfs to the appropriate storage location on GitHub.

For the most part, you can remain blissfully unaware of the details of how this works. For example, if you want to add a new image file to the repository, the process is just like adding any other file to the repository, as long as you have installed git-lfs. In addition, building the documentation does not require any additional steps - again, as long as you have installed git-lfs. But here are some notes that may help you:

.. _obtaining-the-image-files:

### Obtaining the image files

We have configured git-lfs for CTSM so that it does *not* pull down any images by default, since these images are only rarely needed. If you have cloned CTSM or updated an existing clone and want to get the latest version of the image files, you can run:

```shell
git lfs pull --exclude="" --include=""
```

Note that this is done automatically when building the documentation via the `build_docs` command, as described below.

.. _seeing-what-files-are-tracked-by-git-lfs:

### Seeing what files are tracked by git-lfs

To see what file types (extensions) are tracked by git-lfs, run the following from within your CTSM clone:

```shell
git lfs track
```

To see every file currently being managed by git-lfs, run the following from within your CTSM clone:

```shell
git lfs ls-files --exclude="" --include=""
```

.. _adding-a-new-image-file-type:

### Adding a new image file type

If you want to add a new image file type, with an extension not currently being managed by git-lfs, use this process **before** adding the file to the git repository. Note that this assumes that you are in the top-level directory of your CTSM clone:

```shell
git lfs track "*.extension"
git add .gitattributes
git commit -m "Use git-lfs for *.extension files"
```

Then you can check that it worked by running `git lfs track`. Finally, you can add your new file with standard `git add` and `git commit` commands.

.. _procedure-for-building-the-html-documentation:

## Procedure for building the html documentation

.. _recommended-directory-structure-for-these-directions:

### Recommended directory structure for these directions

In order to give concrete examples, we assume a particular directory structure for these directions: We assume you have a directory in your home directory named `ctsm-repos`, and that you clone the main CTSM repository and the ctsm-docs repository inside that directory. So you will have something like:

```
~/ctsm-repos/ctsm
~/ctsm-repos/ctsm-docs
```

It's fine for your `ctsm` directory to be named something more specific, for example, including the branch name as is recommended [on the CTSM GitHub wiki](https://github.com/ESCOMP/CTSM/wiki/Quick-start-to-CTSM-development-with-git). It is also fine to use a different directory organization, but **if you are using the preferred, Docker container-based method, then both your main CTSM clone and the ctsm-docs clone must reside somewhere under your home directory.**

.. _initial-steps-for-building-the-documentation:

### Initial steps for building the documentation

The following procedure assumes you have a clone of CTSM that is checked out at the branch from which you want to build documentation.

Before building the documentation, you need to do the following:

First, run `./bin/git-fleximod update --optional` to update your submodules; this is needed both to get `PTCLM` (which has some files referenced by the documentation build) and to get the `doc-builder` submodule that is used to do the build. **Note the use of the** `--optional` **flag here; this is needed because** `doc-builder` **is an optional submodule.**

Next, if it's not already running, launch the Docker desktop application.

Finally, if you haven't updated the `escomp/base` image in a while, you may want to update to the latest version with:

```shell
docker pull escomp/base
```

.. _building-a-test-version-of-the-documentation-for-your-own-review:

### Building a test version of the documentation for your own review

If you are just building a test version of the documentation for your own review, then you don't need to do anything with the `ctsm-docs` GitHub repository. Instead, you can create a directory according to the [above recommendations for directory structure](#Recommended-directory-structure-for-these-directions).

(However, you may want to follow a procedure similar to the one [documented below for updating the official html documentation](#Updating-the-official-html-documentation): That allows you to do an incremental build rather than building the whole documentation from scratch, which can take about 1/2 hour. You can just avoid committing and pushing when following that procedure, if all you want is to preview the documentation for your own review.)

If you are using the [above recommendations for directory structure](#Recommended-directory-structure-for-these-directions), do the following from the `doc` directory of your ctsm clone:

```shell
mkdir ~/ctsm-repos/ctsm-docs
./build_docs -b ~/ctsm-repos/ctsm-docs -d
```

(where the `-d` flag instructs `build_docs` to use the `escomp/base` Docker container to do the build).

.. _previewing-the-built-documentation:

### Previewing the built documentation

You can then view your changes in a local browser window using a command like one of these:

```shell
open ~/ctsm-repos/ctsm-docs/html/index.html
firefox ~/ctsm-repos/ctsm-docs/html/index.html
konqueror ~/ctsm-repos/ctsm-docs/html/index.html
```

**Note that the version dropdown menu will not work in these local previews.**

.. _updating-the-official-html-documentation:

### Updating the official html documentation

If you want to update the official html documentation, follow the below procedure. (You can also use this procedure if you only intend to build a test version for your own review. This can be helpful because this procedure allows for an incremental build rather than building the entire documentation from scratch.) First clone the ctsm-docs repository to somewhere outside of your main CTSM repository. Here we'll assume that you are using the [above recommendations for directory structure](#Recommended-directory-structure-for-these-directions).

```shell
cd ~/ctsm-repos
git clone https://github.com/ESCOMP/ctsm-docs.git
```

Or, if you already have a copy of ctsm-docs, make sure it is up-to-date as follows:

```shell
cd ~/ctsm-repos/ctsm-docs
git checkout -- .
git pull
```

Next, `cd` to the `doc` directory of your main CTSM clone:

```shell
cd ~/ctsm-repos/ctsm/doc
```

Then, to perform an incremental build (only updating files that need to be updated), run the following. **Note that the following assumes you are building the documentation for the master branch; if you are building the documentation for one of the release branches, replace** `master` **with the appropriate version.**

```shell
./build_docs -b ~/ctsm-repos/ctsm-docs/versions/master -d
```

(where the `-d` flag instructs `build_docs` to use the `escomp/base` Docker container to do the build).

However, if there have been significant changes to the documentation since the last documentation build, and you are intending to push the new build back to GitHub, it's probably best to first clean out the old build in order to remove any no-longer-relevant files. The downside of doing this is that the full documentation rebuild will take about 1/2 hour:

```shell
./build_docs -b ~/ctsm-repos/ctsm-docs/versions/master -d -c
```

(where the `-c` flag instructs `build_docs` to first run `make clean` in the specified directory).

You can preview the documentation as noted [above](#Previewing-the-built-documentation), but now the `index.html` file will be in `~/ctsm-repos/ctsm-docs/versions/master/html/index.html`. **Note that the version dropdown menu appear when you preview the documentation locally, but if you haven't built all the versions, some version links will lead to nonexistent files.**

Then commit and push the rebuilt documentation as follows:

```shell
cd ~/ctsm-repos/ctsm-docs
git add .
git commit -m "YOUR MESSAGE"
git push origin gh-pages
```

**Note if you're using a Windows machine: You may see a lot of warnings like** `warning: LF will be replaced by CRLF in path/to/some/file`. **These seem safe to ignore in this case.**

Within a few minutes, the official documentation at [https://escomp.github.io/ctsm-docs/](https://escomp.github.io/ctsm-docs/) will be updated automatically. Note that you can also fork the ctsm-docs repository and push the gh-pages branch to your fork to allow others to preview the documentation before you overwrite the official documentation.

.. _dealing-with-errors:

### Dealing with errors

.. _inputoutput-error:

#### Input/output error

Some people have reported an error like this:

```shell
Exception occurred:
  File "/usr/lib64/python3.6/shutil.py", line 205, in copystat
    follow_symlinks=follow)
OSError: [Errno 5] Input/output error
The full traceback has been saved in /tmp/sphinx-err-hhpnw_mt.log, if you want to report the issue to the developers.
```

Rerunning the build command will restart it where it left off, and our experience is that that will let you complete the build - though you may need to rerun the build command a few times before it finally finishes successfully.

.. _adding-a-new-version:

### Adding a new version

If you want to add a new version of the built documentation, so that it appears in the dropdown menu, follow this process:

- In the ctsm-docs repository, create a new directory under `versions` with the name of the new version. Our convention is that this name should be the same as the name of the branch in the main CTSM repository that contains this version of the documentation source. 
- In the ctsm-docs repository, add a line in the file `versions/versions.json`. Note that each line is a colon-delimited mapping; the name before the colon is the directory name (i.e., the name of the directory you just created, which is the same as the branch name in the CTSM repository), and the name after the colon is the version name that you want to appear in the dropdown menu. Make sure there is a comma at the end of every line except for the last in this file.
- In the main CTSM repository, check out the branch from which you want to build the documentation. 
- In the main CTSM repository, edit the file `doc/source/conf.py`: edit the `version` and `release` variables (these should probably be the same as the version name that will appear in the dropdown menu) and anything else that needs to be changed for this release (`copyright`, etc.) 
- Build the documentation as described above, being sure to specify your new directory for `BUILDDIR`.

.. _building-a-pdf-of-the-tech-note:

## Building a pdf of the tech note

To build a pdf of the tech note (note that we currently do not support building a pdf of the user's guide), do the following from the `doc` directory of your CTSM clone:

```shell
mkdir ~/ctsm-repos/ctsm-docs-pdf
./build_docs -b ~/ctsm-repos/ctsm-docs-pdf -d -t latexpdf
```

The pdf will appear at `~/ctsm-repos/ctsm-docs-pdf/latex/clmdoc.pdf`.

.. _resources-for-learning-markup-with-restructuredtext-and-using-sphinx:

## Resources for learning markup with reStructuredText and using Sphinx.

* [reStructuredText Primer](http://www.sphinx-doc.org/en/stable/rest.html)
* [ReST Syntax](https://wiki.typo3.org/ReST_Syntax)
* [Sphinx (including how to get Sphinx)](http://www.sphinx-doc.org/en/stable/)
* [reStructured syntax](http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html#tables)

.. _appendix-other-build-methods:

## Appendix: Other build methods

.. _running-build_docs-from-all-of-these-methods:

### Running build_docs from all of these methods

With all of these alternative methods, you can still use `build_docs` commands similar to those given elsewhere in these instructions, **but without the** `-d` **argument**.

.. _relying-more-heavily-on-the-docker-container:

### Relying more heavily on the Docker container

The above instructions still require you to have recent versions of python, git and git-lfs. Alternatively, you can use a method that relies more heavily on the Docker container, avoiding the need for other software to be installed. **This is a very reasonable alternative method, which some people may prefer - especially if you are already comfortable using Docker containers.**

To start, [follow the instructions for installing Docker and downloading the required container](#Install-Docker-and-download-the-required-container). Then use the instructions below.

.. _launching-the-docker-image:

#### Launching the Docker image

Assuming that both your main CTSM repository and ctsm-docs repository are cloned somewhere within your home directory, you can launch the Docker image like this:

```shell
docker run -i -t --rm -v ${HOME}:/home/user/mounted_home escomp/base
```

This tells docker to *run* the container (`docker run`), interactively (`-i`) with a terminal (`-t`), to clean up the container after running (`--rm`), and most importantly, to map your local `${HOME}` directory into the `/home/user/mounted_home` directory inside the container (`-v ${HOME}:/home/user/mounted_home`).

After running this command, you will see a command prompt, but now you will be in the Docker environment. This means:
- You will be able to see directories and files contained in your local `${HOME}` directory, but nothing else. `${HOME}/mounted_home` in the Docker environment refers to `${HOME}` in your native system environment. Changed or added files here will be visible in your native system, in the appropriate subdirectory of your `${HOME}` directory. This will become important later in these directions, since you will view the generated documentation from your native environment.
- Any commands you run from this environment will use the commands bundled with the Docker image, *not* the commands in your native system environment.

When you are done using this Docker image, simply type `exit` from the command prompt.

For far more information on running Docker containers, you can look at [Docker's documentation](https://docs.docker.com/engine/reference/run/).

.. _building-the-documentation:

#### Building the documentation

You can use `build_docs` commands like those given elsewhere in these instructions, **but without the** `-d` **argument** (the `-d` argument launches a new Docker image, but in this case, you are already inside a Docker image, so this is unnecessary). Keep in mind that you will need to specify the path to the documentation build directory according to the Docker file system, not your native file system.

.. _viewing-the-built-documentation:

#### Viewing the built documentation

You can view the built documentation [as described above](#Previewing-the-built-documentation). However, in order to find the `index.html` file, keep in mind the mapping between directories in Docker's environment and those in your native environment.

.. _committing-to-git-repositories:

#### Committing to git repositories

Any local git configuration information is *not* picked up by the Docker image. Therefore, if you want to make any git commits from the Docker terminal session, you will need to do one of the following:
- After doing the `docker run` command as above, copy `/home/user/mounted_home/.gitconfig` to `/home/user/.gitconfig`.
- Or run the `git config` command to set your name and email, as described [here](https://github.com/ESCOMP/CTSM/wiki/Recommended-git-setup#required-settings)

.. _installing-required-software-on-your-desktoplaptop-instructions-for-mac:

### Installing required software on your desktop/laptop (instructions for Mac)

If you choose not to use the Docker-based method described above, you can instead install Sphinx and some related packages locally. **We do not recommend this method, because it can be more challenging to set up and maintain than the method described above. In addition, the instructions here are not regularly tested, so may no longer work. Finally, the use of this method can mean that your Sphinx version differs from the one used to generate the official documentation, which can potentially cause problems.** However, if you cannot use the Docker method for some reason, then you can try this alternative.

.. _prerequisites--assumptions:

#### Prerequisites / assumptions

The following documentation assumes that you have the following available on your machine:

- A recent version of git, ideally configured according to the recommendations in <https://github.com/ESCOMP/CTSM/wiki/Recommended-git-setup#one-time-git-configuration>.

- A recent version of python (python3.5 or later)

- If using a Mac, the [Homebrew package manager](https://brew.sh/)

.. _installing-sphinx-and-the-necessary-sphinx-theme:

#### Installing Sphinx and the necessary Sphinx theme

Sphinx and the necessary Sphinx theme are python packages that can be installed with `pip install`:

```shell
pip install sphinx
pip install git+https://github.com/esmci/sphinx_rtd_theme.git@version-dropdown-with-fixes
```

For more details on installing Sphinx, see the [Sphinx installation docs](http://www.sphinx-doc.org/en/stable/install.html).

*You may then need to update your PATH environment variable to include the path to the* `sphinx-build` *script.*

.. _installing-latexmk:

#### Installing latexmk

`latexmk` is needed both for building the PDF and for creating the equations in the html documentation. On a Mac using homebrew, this can be installed with:

```shell
brew cask install mactex
```

Note that you will need to open a new terminal window for the latexmk tool to be added to your path (it is in `/Library/TeX/texbin/latexmk`, which should have been added to your path by the above `brew` command).

.. _optionally-installing-components-needed-to-build-the-pdf:

#### Optionally installing components needed to build the PDF

If you want to build the PDF (not just the html web pages), you *may* also need rst2pdf. This can be installed with:

```shell
pip install rst2pdf
```

.. _building-the-documentation-on-cheyenne:

### Building the documentation on cheyenne

Finally, although we haven't figured out a way to view the built documentation on cheyenne, you can test the build there, and then transfer the files to your local machine for viewing.

To do so, first ensure that you are using the `git` and `python` modules on cheyenne rather than the default system versions, by doing `module load git` and `module load python`. (Note that the default git module has `git lfs` bundled with it.)

Then, the first time you build the documentation, do this one-time setup:

```shell
pip install --user sphinx
pip install --user sphinxcontrib-programoutput
pip install --user git+https://github.com/esmci/sphinx_rtd_theme.git@version-dropdown-with-fixes
```

and add the following to your `.bashrc` or similar file if using bash:

```bash
export PATH=/glade/u/home/$USER/.local/bin:$PATH
```

or, if you're using tcsh, add to your `.tcshrc`:

```tcsh
setenv PATH /glade/u/home/$USER/.local/bin:$PATH
```

.. _appendix-editing-tips:

## Appendix: Editing tips
- Please don't add manual line breaks when writing text, as this harms searchability. (Note that it's fine to do this in multi-line `:math:` blocks.) For more information, see [text-linebreaks](https://github.com/ESCOMP/CTSM/issues/2135#issuecomment-1764999337).
- You can write the degree symbol ° with Opt-Shift-8 on Mac or Alt+0176 on Windows. Note that this is different from the [masculine ordinal indicator](https://en.wikipedia.org/wiki/Ordinal_indicator) º (typed with Opt-0 on Mac). This is much cleaner and more searchable than using RestructuredText syntax to write a superscript-o!
- Whenever possible, please give equations meaningful labels. E.g., for [eq. 2.26.2](https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/Crop_Irrigation/CLM50_Tech_Note_Crop_Irrigation.html#equation-25-2), `:label: gdds_for_cfts`) instead of numbers (`:label: 25.2`). Numeric labels become obsolete—as you can see in that example!—whenever new equations and/or sections are added. (Note that the equation numbering in the rendered HTML is automatic.)
- Tables defined with the [`:table:` directive](https://docutils.sourceforge.io/docs/ref/rst/directives.html#table) can be annoying because they're very sensitive to the cells inside them being precisely the right widths, as defined by the first `====` strings. If you don't get the widths right, you'll see "Text in column margin" errors. Instead, define your tables using the [`:list-table:`](https://docutils.sourceforge.io/docs/ref/rst/directives.html#list-table) directive.