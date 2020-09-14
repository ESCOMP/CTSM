#!/usr/bin/env bash

# Run emacs/list on all Fortran source files
format_fortran () {
    echo "Parsing $1 as language Fortran"
    emacs --batch -l ./emacs-fortran-formating-script.lisp \
        -f f90-batch-indent-region $1
}
export -f format_fortran
find ../lilac/ -iregex ".*\.F[0-9]*" -exec bash -c 'format_fortran "$0"' {} \;
find ../tests/ -iregex ".*\.F[0-9]*" -exec bash -c 'format_fortran "$0"' {} \;
