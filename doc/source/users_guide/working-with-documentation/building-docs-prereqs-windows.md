.. _building-docs-prereqs-windows:

# Initial setup: Windows

Note that you may need administrator privileges on your PC (or approval from your IT department) for various steps here.

.. _install-wsl:

## Install Linux subsystem

We don't support building our documentation in the native Windows command-line environment. Thus, you will need to install a little version of Linux inside a virtual machine (VM) to use instead.

1. Download and install Ubuntu from the Microsoft Store.
1. Restart your computer.
1. Open Ubuntu.

If Ubuntu opens in that last step but you see an error, you may need to manually enable Windows Subsystem for Linux (WSL). To do so: Open Control Panel, go to "Programs" > "Programs and Features" > "Turn Windows features on or off". Check the box next to "Windows Subsystem for Linux" and click OK.

.. _windows-docs-ubuntu-utilities:

## Install utilities
Enter the following commands **into your Ubuntu terminal** to install any missing utilities we need (the `which ... ||` should make it so that no installation happens if you already have it):
```shell
# Refresh the list of available software
sudo apt-get update

# make: Part of the docs-building process
which make || sudo apt-get -y install make

# git and git-lfs, needed for getting and contributing to the CTSM code and docs
which git || sudo apt-get -y install git
which git-lfs || sudo apt-get -y install git-lfs

# Chromium: A web browser engine that's the basis for popular browsers like Google
# Chrome and Microsoft Edge
which chromium || sudo apt-get -y install chromium
```

.. _container-or-conda-windows:

## Install container software or Conda environment

We recommend building the software in what's called a container—basically a tiny little operating system with just some apps and utilities needed by the doc-building process. This is nice because, if we change the doc-building process in ways that require new versions of those apps and utilities, that will be completely invisible to you. You won't need to manually do anything to update your setup to work with the new process; it'll just happen automatically.

We recommend using the container software Podman.

1. Install Podman with `sudo apt-get -y install podman`.
1. Set up and start a Podman "virtual machine" with `podman machine init --now`.
1. Test your installation by doing `podman run --rm hello-world`. If it worked, you should see ASCII art of the Podman logo.

You may not be able to install Podman or any other containerization software, so there is an alternative method: a Conda environment.

1. Check whether you already have Conda installed by doing `which conda`. If that doesn't print anything, [install Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux).
1. Follow the instructions for setting up the `ctsm_pylib` Conda environment in Sect. :numref:`using-ctsm-pylib`.


## Set up your permissions
This will make sure that you "own" your home directory in the Ubuntu VM. **In your Ubuntu terminal**, do:
```shell
chown -R $USER:$USER $HOME
```

.. _editing-text-files-wsl:

## Editing text files in an Ubuntu VM
If you prefer using an old-school text editor like `vim`, it's probably already installed, or can be installed with `sudo apt-get -y install EDITOR_NAME`. If you prefer a more user-friendly interface, there are several options.

You may be able to edit files in your Ubuntu VM in the Ubuntu terminal by using the name of the Windows executable. For Notepad, for instance, you would do 
```shell
notepad.exe file_i_want_to_edit.rst
```

If you use [VS Code](https://code.visualstudio.com/), you can install the [WSL VS Code extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl). Then you can open any file or folder in your Ubuntu VM by doing
```shell
code path/to/file-or-folder
```

You can also install a user-friendly text editor in Ubuntu. This may be slower and have unexpected differences in behavior from what you expect from Windows apps, but it does work. For example:
- [gedit](https://gedit-text-editor.org/): `sudo apt-get install -y gedit`
- [Kate](https://kate-editor.org/): `sudo apt-get install -y kate`
- [VS Code](https://code.visualstudio.com/) (if you don't already have it installed on Windows): `sudo snap install code --classic`

You can use all of those to open and edit files, but Kate and VS Code let you open entire folders, which can be convenient. In any case, you'd do `EDITOR_NAME path/to/thing/youre/editing` to open it, where `EDITOR_NAME` is `gedit`, `kate`, or `code`, respectively.
