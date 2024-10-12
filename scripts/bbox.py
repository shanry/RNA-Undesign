#!/usr/bin/env python3
# need pip3 install PyMuPDF (not pip3 install fitz!)

import sys
import fitz

pdf = fitz.open(sys.argv[1])
p = next(pdf.pages()) # first page
bbox = p.bleedbox
print(abs(bbox.x0 - bbox.x1), abs(bbox.y0 - bbox.y1))
