#!/bin/bash
# dependency: `pdfjam pdfcrop`
# usage: `./pdfrotate.sh path.pdf` generates pdf fils for each plotstr in plotstrs.txt

pdf_path=$1
echo "pdf path:" $pdf_path

# Get the basename (file name)
file_name=$(basename "$pdf_path")
echo "File name: $file_name"

# Get the directory path
dir_name=$(dirname "$pdf_path")
echo "Directory path: $dir_name"


# code for automatic rotation without pdfjam
rm -f /tmp/angles
# INPUT_PDF="${pdf_path}"
# INPUT_PDF_ABS_PATH="$(readlink -f $INPUT_PDF)"
for angle in `seq 0 30 179` # every 30 degrees
do
  echo "rotating", $angle
  pdfjam --angle $angle --outfile /tmp/${file_name}-${angle}.pdf $pdf_path > /dev/null
  pdfcrop /tmp/${file_name}-${angle}.pdf > /dev/null
  (echo -ne "$angle\t"; ./scripts/bbox.py /tmp/${file_name}-${angle}-crop.pdf) >> /tmp/angles
done

bestangle=`cat /tmp/angles | sort -nk3 | head -1 | cut -f 1` # smallest height (y)
echo crude best angle, $bestangle

# exit 0

low=$((bestangle-20))
high=$((bestangle+20))

echo trying refined angles from $low to $high

for angle in `seq $low 5 $high` # every 5 degrees
do
  echo "rotating", $angle
  pdfjam --angle $angle --outfile /tmp/${file_name}-${angle}.pdf $pdf_path > /dev/null
  pdfcrop /tmp/${file_name}-${angle}.pdf > /dev/null
  (echo -ne "$angle\t"; ./scripts/bbox.py /tmp/${file_name}-${angle}-crop.pdf) >> /tmp/angles
done

bestangle=`cat /tmp/angles | sort -nk3 | head -1 | cut -f 1` # smallest height (y)
echo refined best angle, $bestangle

mv /tmp/${file_name}-${bestangle}-crop.pdf "outputs/${file_name}"
echo "final output: outputs/${file_name}"
echo "------------------------"
echo