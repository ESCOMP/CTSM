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

## Install utilities
Enter the following commands **into your Ubuntu terminal** to install any missing utilities we need (the `which ... ||` should make it so that no installation happens if you already have it):
```shell
sudo apt-get update

# Podman: The software that will run the container in which the docs are built
which podman || sudo apt-get -y install podman

# make: Part of the docs-building process
which make || sudo apt-get -y install make

# git and git-lfs, needed for getting and contributing to the CTSM code and docs
which git || sudo apt-get -y install git
which git-lfs || sudo apt-get -y install git-lfs
```

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
