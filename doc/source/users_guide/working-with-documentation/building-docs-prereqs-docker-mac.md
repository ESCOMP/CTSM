.. _building-docs-prereqs-docker-mac:

# Initial setup: Mac

Note that you may need administrator privileges on your Mac for the installation steps detailed here.

.. _building-docs-git-tools:

## Python
To test whether you already have the required Python version, open a Terminal window and try the following:
```shell
python3 --version
```

If python3 is already set up, you'll see a version number. If that version is 3.7 or later, you should be ready as far as Python goes; continue to :ref:`docs-git-tools`.

If not, recent versions of macOS should print a messsage saying, "xcode-select: No developer tools were found, requesting install." A dialog box will then pop up that says, "The 'python3' command requires the command line developer tools. Would you like to install the tools now?" Press Install and go through the installation process. This will take a while; once it's done, test by doing ``python3 --version`` again. (You may need to open a new Terminal window.) If the printed version number looks good, continue to :ref:`docs-git-tools`.

..
  The paragraph above was tested 2025-04-25 on a fresh-ish installation of macOS 15.3.2.

If instead `python3` gives "command not found," or the version is less than 3.7, you might need to install Python; continue to :ref:`aliasing-python3-to-python`.

.. _aliasing-python3-to-python:

### Aliasing `python3` to `python`
Try the same command as above, but instead of `python3` just do `python` (no number). If that version is 3.7 or later, you can tell your Mac that when you say `python3` you want it to use `python`:
```bash
echo alias python3="$(which python)" >> ~/.bashrc
echo alias python3="$(which python)" >> ~/.zshrc
```

This will make it so that bash scripts, like what we use to build our docs, know what to do for `python3`. `python3` will also be available in new Terminal sessions if your shell is `zsh` (the default since macOS 10.15) or `bash`.

If you were able to do this, you can continue to :ref:`docs-git-tools`.

### Conda
If your `python` doesn't exist or is too old, we suggest using Python via Conda. First, check whether you already have Conda installed:
```shell
conda env list
```

If that shows you something like
```
# conda environments:
#
base           /Users/you/...
another_env    /Users/you/.../...
...
```

instead of the "command not found" error, then you do have conda installed! (Note that the second column doesn't really matter.) Try this to check the Python version in the `base` Conda environment:
```shell
conda run -n base python3 --version
```

Repeat with all your Conda environments as needed until you find one that's Python 3.7 or later. Let's say your `ENVNAME` environment works. In that case, just make sure to do `conda activate ENVNAME` before running the commands in the documentation-building instructions.

If you don't have Conda yet, go on to the next section. Otherwise, continue to :ref:`docs-git-tools`.

.. _installing-conda-for-docs:

### Installing Conda
We suggest installing Conda, if needed, via Miniforge:

1. [Download Miniforge](https://conda-forge.org/download/) and install it. (:ref:`what-kind-of-mac-chip`) You can also [install Miniforge via Homebrew](https://formulae.brew.sh/cask/miniforge#default), if you already have that installed.
2. Activate Conda permanently in your shell by opening a new Terminal window and doing `conda init "$(basename $SHELL)"`.

You should now have `conda` and an up-to-date version of `python3` available, although will need to open another new Terminal window for it to work.

## Additional requirements
The remaining software you need to install depends on which documentation-building method(s) you plan to use.

### Docker (recommended)
Install Docker Desktop through the [download page on Docker's website](https://docs.docker.com/desktop/setup/install/mac-install/). (:ref:`what-kind-of-mac-chip`)

### Non-Docker method

Install Conda, if needed (see :ref:`installing-conda-for-docs`). Then follow the instructions for setting up the `ctsm_pylib` Conda environment in Sect. :numref:`using-ctsm-pylib`.

.. _docs-git-tools:

## Git tools
Note: Do this section after handling Python, because the Python installation process might bring the Git tools with it.

To test whether you have the required Git tools already, open a Terminal window and try the following:
```shell
git --version
git-lfs --version
```

If either of those fail with "command not found," you'll need to install them. The recommended way is with [Homebrew](https://brew.sh/).

1. [Install Homebrew](https://brew.sh/), if needed. (Check using `brew --version`.) Make sure to follow the instructions for adding Homebrew to your path. (:ref:`what-kind-of-mac-chip`)
2. Use Homebrew to [install Git](https://formulae.brew.sh/formula/git#default), if needed.
3. Use Homebrew to [install Git LFS](https://formulae.brew.sh/formula/git-lfs#default), if needed.

## Frequently-asked questions

.. _what-kind-of-mac-chip:

### What kind of chip does my Mac have?
For certain steps in this installation process, you may need to know whether your Mac has an Intel (`x86_64`) or an Apple Silicon (`arm64`) chip. If you don't know, visit Apple's [Mac computers with Apple silicon](https://support.apple.com/en-us/116943) page for instructions.