.. _building-docs-prereqs-docker-windows:

# Initial setup: Windows

Note that you may need administrator privileges on your PC (or approval from your IT department) for various steps here.

## Linux subsystem

We don't support building our documentation in the native Windows command-line environment. Thus, you will need to install a little version of Linux to use instead.

The first step is to install Windows Subsystem for Linux (WSL), if you haven't already:

1. Download and install it from the Microsoft Store.
1. Restart your computer.

Next, download and install Ubuntu (a version of Linux) via the Microsoft Store.

Finally, try opening Ubuntu. If you run into an error, you may need to manually enable WSL. To do so: Open Control Panel, go to "Programs and Features," and then "Turn Windows features on or off." Check the box next to "Windows Subsystem for Linux" and click OK.

## Docker Desktop
1. Install Docker Desktop through the Microsoft Store or the [download page on Docker's website](https://docs.docker.com/desktop/setup/install/windows-install/).
1. Open Docker Desktop. In Settings > Resources > WSL, Enable Docker to work under Ubuntu.
