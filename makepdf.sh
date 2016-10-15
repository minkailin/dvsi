#!/bin/bash
latex paper.tex
bibtex paper
latex paper
latex paper
dvips -o paper.ps paper.dvi
ps2pdf -sPAPERSIZE=a4 paper.ps
echo 'done'
