#!/bin/bash
# dependency: `pip install pdfCropMargins`
# usage: `cat plotstrs.txt | xargs -L 1 ./draw_motif.sh` generates pdf fils for each plotstr in plotstrs.txt

# Check if the environment variable is not set
if [ -z "${VIENNA}" ]; then
    # Set the variable if it's not set
    VIENNA=~/biology/ViennaRNA-2.5.1
fi
echo "\$VIENNA: $VIENNA"
mode=0
# if [ -z "${PLMD}" ]; then
#     mode=0
# else
#     mode=$PLMD
# fi
echo "mode:" $mode

# read -r line
line=$1
echo $line

# drawfunction=${2:-drawpoints}
echo "drawfunction:" $drawfunction

ID=$(echo "$line" | cut -d',' -f1) #  ID
prestring=$(echo "$line" | cut -d',' -f3) # annotation
# poststring=$(echo "$line" | cut -d',' -f4) # annotation

# Check if the file already exists
# order=0
# while [ -e "${id}${order}_mode${mode}.pdf" ]; do
#     echo "exist: ${id}${order}_mode${mode}.pdf"
#     ((order++))
# done
# id="${id}${order}_mode${mode}"
# fi
id="${ID}_automode${mode}"
echo "id:" $id

struct=$(echo "$line" | cut -d',' -f2) # structure
len=${#struct}
echo length $len
seq=`for i in \`seq $len\`; do echo -ne o; done` # o....o for seq
echo "seq   " $seq
echo struct $struct

prestring=$(echo "$line" | cut -d',' -f3) # annotation
poststring=$(echo "$line" | cut -d',' -f4) # annotation

echo "prestring " $prestring
echo "poststring " $poststring

# Initialize an empty array to store the transformed second indices
transformed_indices=()

# Use awk to extract the first and second indices and subtract 1 from the second one
while read -r first second; do
    transformed_first=$((first - 1))    # Subtract 1 from the first index
    transformed_second=$((second - 1))  # Subtract 1 from the second index
    transformed_indices+=("$transformed_first" "$transformed_second")
done < <(echo "$poststring" | awk '{for(i=1; i<=NF; i+=5) print $i, $(i+1)}')

# Format the output in PostScript style
skip_string="/my_list [${transformed_indices[*]}] def"
echo "skip list:" $skip_string

# exit
if [ -z "${PATH_FASTMOTIF}" ]; then
    PATH_BASE=`pwd`
else
    echo "PATH_FASTMOTIF:" $PATH_FASTMOTIF
    PATH_BASE=$PATH_FASTMOTIF
fi
echo "PATH_BASE:" $PATH_BASE

# # -t 0 means layout mode 0 (default 1)
echo -ne ">$id\n$seq\n$struct" | $VIENNA/bin/RNAplot -t $mode --pre "$prestring " --post "$poststring" # the final space is important to keep "" #"$span GREEN Fomark"

# Call the Python script and capture its output
result_cross=$(python ${PATH_BASE}/scripts/coord.py ${id}_ss.ps)
echo "result_cross:" $result_cross

if [ "$result_cross" -eq 1 ]; then
    echo "cross"
    # mv ${id}_ss.ps cross.ps
    mode=4
    id="${ID}_automode${mode}"
    echo -ne ">$id\n$seq\n$struct" | $VIENNA/bin/RNAplot -t $mode --pre "$prestring " --post "$poststring" # the final space is important to keep "" #"$span GREEN Fomark"
    # python ${PATH_BASE}/scripts/replace.py ${id}_ss.ps cross.ps
    # cp cross.ps ${id}_ss.ps
fi

sed -i 's/fsize setlinewidth/8 setlinewidth/' ${id}_ss.ps # change line width
sed -i '/\/colorpair/,/grestore/{s/hsb/1.0\n  sethsbcolor\n  3 pop/}' ${id}_ss.ps # change color
sed -i 's/0.667 0.5 colorpair/0.583 1.0 colorpair/g' ${id}_ss.ps # change color
sed -i 's/0.70 0.5 colorpair/0.583 1.0 colorpair/g' ${id}_ss.ps # change color
sed -i 's/0.1667 1.0 colorpair/0.1083 1.0 colorpair/g' ${id}_ss.ps # change color
cp ${id}_ss.ps ${id}.ps


sed "/^init$/ {
    r ${PATH_BASE}/scripts/insert.ps
}
s/^drawoutline$/drawarrows\ndrawpoints/
# s/^drawoutline$/drawpoints\ndrawarrows/
# s/^drawoutline$/drawarrows/
s/^drawbases$//
s/BoundingBox: 0 0/BoundingBox: -10 -10/
## disable drawpairs
# s/^drawpairs$//

" ${id}.ps > ${id}_ss.ps

sed -i "s|/my_list \[0\] def|$skip_string|" ${id}_ss.ps # skip list

rm ${id}.ps

substring="1 1 10 WHITE"

target_str="0 1 coor_len 1 sub"
new_str="1 1 coor_len 2 sub"

# Check if the file contains substring_a
if [[ "$poststring" == *"$substring"* ]]; then
    # If substring_a is found, replace target_str with new_str in the file
    sed -i "s/$target_str/$new_str/g" "${id}_ss.ps"
    sed -i "s/i 0 eq/i -2 eq/g" "${id}_ss.ps"
    sed -i "s/i coor length 1 sub eq/i -2 eq/g" "${id}_ss.ps"
    echo "Replaced '$target_str' with '$new_str' in the file."
else
    echo "Substring '$substring' not found in the file."
fi

ps2pdf -dEPSCrop ${id}_ss.ps > /dev/null # bounding box
# mv ${id}_ss.pdf ${id}.pdf
# crop margin automatically; alternative way: pdfcrop ${id}_ss.pdf
pdfcropmargins -v -u -s ${id}_ss.pdf -o ${id}_ss-crop.pdf > /dev/null # pip install pdfCropMargins
mv ${id}_ss-crop.pdf ${id}.pdf # final output
rm ${id}_ss.p* # remove temp files

mkdir -p outputs
mv ${id}.pdf outputs/ # move to outputs folder

echo "outputs/${id}.pdf"
exit
