```shell
cd doc
./build_docs -b _build -c -d
```

This runs a complicated series of scripts and software that culminate in something called Sphinx converting the .rst and .md files into HTML webpages, all within the (new, if needed) `_build` directory.

The `-c` means "do a clean build." If you leave it off, Sphinx will only rebuild files it thinks have changed since the last time you built the docs, which is faster but can lead to problems. If you get unexpected results or errors without `-c`, rerunning with `-c` is the first troubleshooting step.

The `-d` means "run using the container." If you're not using the container and are instead using the `ctsm_pylib` Conda environment **(not recommended)**, leave off `-d`, and make sure you've activated `ctsm_pylib` before running the above command.

(Do `./build_docs --help` for more information and options.)
