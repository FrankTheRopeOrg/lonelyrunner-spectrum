#!/bin/sh
# usage: ./set_doi.sh 10.5281/zenodo.1234567 [path/to/lonelyrunner.tex]
# Replaces the placeholder DOI in CITATION.cff, README.md and (optionally) the paper.
[ -z "$1" ] && { echo "usage: $0 DOI [tex-file]"; exit 1; }
DOI="$1"
sed -i "s#10.5281/zenodo.XXXXXXX#$DOI#g" CITATION.cff README.md
[ -n "$2" ] && sed -i "s#10.5281/zenodo.XXXXXXX#$DOI#g" "$2"
grep -n "$DOI" CITATION.cff README.md ${2:+"$2"}
