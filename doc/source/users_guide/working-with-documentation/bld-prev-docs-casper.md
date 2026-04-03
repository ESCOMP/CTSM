(bld-prev-docs-casper)=

# Building and previewing the documentation on Casper (RECOMMENDED)

```{contents}
:depth: 1
:local:
```

## Initial Casper setup

None! This is why Casper is the recommended method for building and previewing the docs.

## Building docs on Casper
Casper uses the Podman software for running containers like the one we recommend for building the CTSM documentation. Make sure it's enabled before building: `module load podman`. Then do:

```{include} embed-build-cmd.md
```

See the "Container software or Conda environment" sections for {ref}`Mac <container-or-conda-mac>` or {ref}`Windows <container-or-conda-windows>` for more information on these two methods.

## Previewing docs on Casper

There are a few different ways to do this.

### Previewing docs on Casper with OnDemand

The simplest way to preview your built documentation on Casper is to use [NCAR's OnDemand service](https://ondemand.hpc.ucar.edu/pun/sys/dashboard). After you [open a new Casper Login VNC Desktop session](https://ondemand.hpc.ucar.edu/pun/sys/dashboard/batch_connect/sys/login_desktop_ncar/session_contexts/new), wait for it to start, then click the "Launch Casper Login VNC Desktop" button. This will open a Linux desktop in your browser. Click on the "Home" icon and navigate to your CTSM checkout, then the `doc/_build/html` directory, then open the `index.html` file. That should open it in the Linux desktop's Firefox, at which point you can browse around as usual.

If you rebuild the documentation but don't see your changes updated in the webpage, and the reload button doesn't work, you may need to navigate to a different page and then come back.

### Previewing docs on Casper via SSH tunnel

If the OnDemand method doesn't work for whatever reason, open a terminal on your local machine and do this:
```shell
ssh YOUR_USERNAME@casper.hpc.ucar.edu echo $((10000 + $(id -u) % 50000))
```

This will print an integer that we will call `YOUR_PORT`.

Then open a new SSH connection to the server like so, replacing `YOUR_PORT` with the integer you got above:
```shell
ssh -L YOUR_PORT:localhost:YOUR_PORT YOUR_USERNAME@casper.hpc.ucar.edu
```

Once that's connected, we're going to spin up a web server. It's best to do this on a compute node rather than a login node. Use the `qinteractive` command to open an interactive session on a Casper compute node. (Find more info about `qinteractive` [here](https://ncar-hpc-docs.readthedocs.io/en/latest/pbs/#qinteractive).)

Once your interactive compute session is open, start the web server in it like so (again replacing `YOUR_PORT`):
```shell
cd /path/to/your/ctsm/repo
cd doc/_build/html
python3 -m http.server YOUR_PORT
```

Now you're ready to view your documentation! Just open a web browser on your computer and navigate to `http://localhost:YOUR_PORT`. You should see a rendered version of the documentation that you can browse as usual.

### Fallback preview method

If neither of the above work, you can download the `doc/_build/html` directory to your computer (e.g., with `scp` or `rsync`) and open `html/index.html` in a web browser.

### A note about previewing the docs

```{include} embed-preview-menu.md
```
