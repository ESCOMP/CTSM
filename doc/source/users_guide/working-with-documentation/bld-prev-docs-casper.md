.. _bld-prev-docs-casper:

# Building and previewing the documentation on Casper (RECOMMENDED)

.. contents::
   :depth: 1
   :local:

## Initial Casper setup

None! This is why Casper is the recommended method for building and previewing the docs.

## Building docs on Casper
Casper uses the Podman software for running containers like the one we recommend for building the CTSM documentation. Make sure it's enabled before building: `module load podman`. Then do:

.. mdinclude:: embed-build-cmd.md

See the "Container software or Conda environment" sections for :ref:`Mac <container-or-conda-mac>` or :ref:`Windows <container-or-conda-windows>` for more information on these two methods.

## Previewing docs on Casper

Open a terminal on your local machine and do this:
```shell
# user@server e.g. samrabin@casper.hpc.ucar.edu
ssh user@server echo $((10000 + $(id -u) % 50000))
```

This will print an integer that we will call `YOUR_PORT`.

Then open a new SSH connection to the server like so, replacing `YOUR_PORT` with the integer you got above:
```shell
ssh -L YOUR_PORT:localhost:YOUR_PORT user@server  # e.g., samrabin@casper.hpc.ucar.edu
```

Once that's connected, we're going to spin up a web server. It's best to do this on a compute node rather than a login node. For Casper, use the `qinteractive` command to open an interactive session on a compute node. (Find more info about `qinteractive` [here](https://ncar-hpc-docs.readthedocs.io/en/latest/pbs/#qinteractive).)

Once your interactive compute session is open, start the web server in it like so (again replacing `YOUR_PORT`):
```shell
cd /path/to/your/ctsm/repo
cd doc/_build/html
python3 -m http.server YOUR_PORT
```

Now you're ready to view your documentation! Just open a web browser on your computer and navigate to `http://localhost:YOUR_PORT`. You should see a rendered version of the documentation that you can browse as usual.

If that doesn't work, you can download the `doc/_build/html` directory to your computer (e.g., with `scp` or `rsync`) and open `html/index.html` in a web browser.

.. mdinclude:: embed-preview-menu.md
