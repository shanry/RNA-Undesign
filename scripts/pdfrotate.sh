#!/bin/bash
# dependency: `pdfjam pdfcrop`
# usage: `./pdfrotate.sh path.pdf` generates pdf fils for each plotstr in plotstrs.txt

pdf_path=$1
echo "pdf path:" $pdf_path

flag_angle=0
if [ $# -gt 1 ]; then
  angle=$2
  flag_angle=1
else
  flag_angle=0
fi
echo "angle:" $angle

# Get the basename (file name)
file_name=$(basename "$pdf_path")
echo "File name: $file_name"

# Get the directory path
dir_name=$(dirname "$pdf_path")
echo "Directory path: $dir_name"

# Rotate the pdf according to the user-defined angle
if [ $flag_angle -eq 1 ]; then
  pdfjam --angle $angle --outfile ${file_name}-${angle}.pdf $pdf_path > /dev/null
  echo "------------------------"
  echo "final output: ${file_name}-${angle}.pdf"
  exit 0
fi

exit 0

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