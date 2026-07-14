#!/bin/bash

set -e

TEX_FILE="tex/SpectralO33.tex"
OUTPUT_DIR="out"
MAIN_NAME="SpectralO33"

mkdir -p "$OUTPUT_DIR"
export TEXINPUTS=".:./tex:${TEXINPUTS}"

pdflatex -file-line-error -interaction=nonstopmode -synctex=1 \
  -output-directory="$OUTPUT_DIR" -recorder "$TEX_FILE"

cd "$OUTPUT_DIR"
BSTINPUTS="../tex:${BSTINPUTS}" BIBINPUTS="../tex:${BIBINPUTS}" bibtex "$MAIN_NAME"
cd ..

pdflatex -file-line-error -interaction=nonstopmode -synctex=1 \
  -output-directory="$OUTPUT_DIR" -recorder "$TEX_FILE"
pdflatex -file-line-error -interaction=nonstopmode -synctex=1 \
  -output-directory="$OUTPUT_DIR" -recorder "$TEX_FILE"

test -f "$OUTPUT_DIR/$MAIN_NAME.pdf"
echo "Compiled $OUTPUT_DIR/$MAIN_NAME.pdf"
