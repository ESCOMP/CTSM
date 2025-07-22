.. _building-docs-prereqs-windows:

# Initial setup: Windows

Note that you may need administrator privileges on your PC (or approval from your IT department) for various steps here.

.. _install-wsl:

## Install Linux subsystem

We don't support building our documentation in the native Windows command-line environment. Thus, you will need to install a little version of Linux inside a virtual machine (VM) to use instead. The process for doing this varies depending on how tightly the installation process is controlled on your computer.

### NCAR computers

Please follow the [Windows Subsystem for Linux (WSL) setup instructions](https://wiki.ucar.edu/pages/viewpage.action?pageId=514032264&spaceKey=CONFIGMGMT&title=Setup) on the UCAR Wiki. In the step about installing a Linux distribution, choose Ubuntu.

Feel free to peruse the [overall WSL documentation](https://wiki.ucar.edu/spaces/CONFIGMGMT/pages/514032242/Windows+Subsystem+for+Linux) on and linked from the UCAR Wiki for additional information.

### Non-NCAR computers

If your computer is managed by an organization other than NCAR, please check with your IT department or equivalent for instructions on installing Windows Subsystem for Linux (WSL) and Ubuntu. Otherwise, follow these instructions:

1. Download and install Ubuntu from the Microsoft Store.
1. Restart your computer.
1. Open Ubuntu.

If Ubuntu opens in that last step but you see an error, you may need to manually enable Windows Subsystem for Linux (WSL). To do so: Open Control Panel, go to "Programs" > "Programs and Features" > "Turn Windows features on or off". Check the box next to "Windows Subsystem for Linux" and click OK.

Once Ubuntu is working and open, you'll be asked to create a new UNIX username and password. This doesn't have to match your Windows username and password, but do make sure to save this information somewhere secure.

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

# WSL utilities, which will give us the wslview command for opening HTML pages in a Windows browser
which wslview || sudo apt-get -y install wslu
```

.. _container-or-conda-windows:

## Install container software or Conda environment

We recommend building the software in what's called a container—basically a tiny little operating system with just some apps and utilities needed by the doc-building process. This is nice because, if we change the doc-building process in ways that require new versions of those apps and utilities, that will be completely invisible to you. You won't need to manually do anything to update your setup to work with the new process; it'll just happen automatically.

For builds in WSL (Ubuntu), we recommend using the container software Docker. You can install it in Ubuntu like so:

```shell
# If needed, download and run the Docker installation script.
# Ignore the message saying "We recommend using Docker Desktop for Windows."
# The script will make you wait 20 seconds to make sure this is want you want,
# and then it should continue automatically.
which docker || curl -fsSL https://get.docker.com -o get-docker.sh
which docker || sudo sh ./get-docker.sh

# Set up the docker "group," if needed, and add your username to it.
sudo groupadd docker  # Create docker group if it doesn't exist
sudo usermod -aG docker $USER  # Add your user to the docker group
newgrp docker  # Apply the new group membership (avoids needing to log out and back in)

# Make sure it worked: This should print a "Hello from Docker!" message
docker run hello-world
```

You may not be able to install Docker or any other containerization software, so there is an alternative method: a Conda environment.

1. Check whether you already have Conda installed by doing `which conda`. If that doesn't print anything, [install Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux).
1. Follow the instructions for setting up the `ctsm_pylib` Conda environment in Sect. :numref:`using-ctsm-pylib`.

.. _editing-text-files-wsl:

## Editing documentation files
If you prefer using an old-school text editor like `vim`, it's probably already installed in your Ubuntu VM, or can be installed with `sudo apt-get -y install EDITOR_NAME`. If you prefer a more user-friendly interface, there are several options. Note that **all commands in this section are to be run in your Ubuntu VM, not a Windows terminal**.

### In a Windows app (recommended)
If you installed `wslview` in the instructions above, you can edit files by doing 
```shell
wslview path/to/file_i_want_to_edit.rst
```
If not, you can do
```shell
explorer.exe $(wslpath -w path/to/file_i_want_to_edit.rst)
```
These both do the same thing, but the `wslview` method is simpler. Either way, at least the first time you do this, it will open a window asking which app you'd like to open the file in. Choose whatever you're most comfortable with. At the bottom of the window, you can then choose whether you always want to open HTML files using the selected app or just this once.

You may also be able to open files in Windows apps by using the name of the Windows executable. For Notepad, for instance, you would do 
```shell
notepad.exe $(wslpath -w path/to/file_i_want_to_edit.rst)
```

If you use [VS Code](https://code.visualstudio.com/), you can install the [WSL VS Code extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl). Then (after closing and re-opening Ubuntu) you can open any documentation file **or folder** by doing
```shell
code path/to/file-or-folder
```

### In an Ubuntu app (not recommended)

You can also install a user-friendly text editor in Ubuntu. This may be slower and have unexpected differences in behavior from what you expect from Windows apps, but it does work. For example:
- [gedit](https://gedit-text-editor.org/): `sudo apt-get install -y gedit`
- [Kate](https://kate-editor.org/): `sudo apt-get install -y kate`
- [VS Code](https://code.visualstudio.com/) (if you don't already have it installed on Windows): `sudo snap install code --classic`

You can use all of those to open and edit files, but Kate and VS Code let you open entire folders, which can be convenient. In any case, you'd do `EDITOR_NAME path/to/thing/youre/editing` to open it, where `EDITOR_NAME` is `gedit`, `kate`, or `code`, respectively.

## Troubleshooting

### "Permission denied" error

If you get this error, it may be a result of opening Ubuntu as an administrator (e.g., by right-clicking on its icon and choosing "Run as administrator.") Try not doing that, although this will result in you needing to get a new copy of CTSM to work in. 

If that's not feasible or doesn't solve the problem, you may need to remind Linux that you do actually own your files. **In your Ubuntu terminal**, do:
```shell
chown -R $USER:$USER $HOME
```

If that also gives a permission error, you may need to put `sudo` at the start of the command.

### "The host 'wsl$' was not found in the list of allowed hosts"

You may see this warning in a dialog box after trying to open a file with `wslview`, `explorer.exe`, or something else. Check "Permanently allow host 'wsl$'" and then press "Allow".
