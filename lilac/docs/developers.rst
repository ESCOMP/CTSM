Developers Guide to Using LILAC
===============================

Building LILAC
--------------

LILAC can be build using CMake::

    $ cd /lilac/build && cmake ..
    $ make

For development and testing purposes, LILAC can also be built using a
`docker-compose` script::

    $ docker-compose build
    $ docker-compose run

Testing LILAC
-------------

LILAC includes a full test suite including unit tests and a number of
simplified coupled model tests. To run these tests, simply run the `ctest`
command::

    $ ctest

Note, if you are using the docker-compose development, the `docker-compose run`
command will build and run LILAC automatically.
