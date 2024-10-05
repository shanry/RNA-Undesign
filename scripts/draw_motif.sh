#!/bin/bash
# dependency: `pip install pdfCropMargins`
# usage: `cat plotstrs.txt | xargs -L 1 ./plot_motif.sh` generates pdf fils for each plotstr in plotstrs.txt

# Check if the environment variable is not set
if [ -z "${VIENNA}" ]; then
    # Set the variable if it's not set
    VIENNA=~/biology/ViennaRNA-2.5.1
fi
echo "\$VIENNA: $VIENNA"

mode=0
echo "mode:" $mode

# read -r line
line=$1
echo $line

id=$(echo "$line" | cut -d',' -f1) #  ID
prestring=$(echo "$line" | cut -d',' -f3) # annotation
# poststring=$(echo "$line" | cut -d',' -f4) # annotation

# Check if the file already exists
order=0
while [ -e "${id}${order}_mode${mode}.pdf" ]; do
    echo "exist: ${id}${order}_mode${mode}.pdf"
    ((order++))
done
id="${id}${order}_mode${mode}"
# fi
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

# exit

# paircount=$(echo "$line" | cut -d',' -f4) # eterna ID
# echo "paircount:" $paircount
# # prestring="$first $second GREEN Fomark"

# # for ((c=0; c<paircount; c++)) do
# # 	i=$(echo "$line" | cut -d',' -f$((4+2*c)))
# # 	j=$(echo "$line" | cut -d',' -f$((4+2*c+1)))
# # 	prestring="$prestring $(($i+1)) $(($j+1)) 0 0.6 colorpair"
# # 	echo "i:" $i  "j:" $j
# # 	echo $prestring
# # done
# # exit 0
# if [ $paircount -eq 1 ]; then
# 	echo "one pair"
# 	i=$(echo "$line" | cut -d',' -f5)
# 	j=$(echo "$line" | cut -d',' -f6)
# 	echo "i:" $i  "j:" $j
# 	prestring="$(($i+1)) $(($j+1)) 0 5 8 Fomark"
# elif [ $paircount -eq 2 ]; then
# 	echo "two pairs"
# 	i=$(echo "$line" | cut -d',' -f5)
# 	j=$(echo "$line" | cut -d',' -f6)
# 	k=$(echo "$line" | cut -d',' -f7)
# 	l=$(echo "$line" | cut -d',' -f8)
# 	echo "i:" $i  "j:" $j  "k:" $k  "l:" $l
# 	prestring="$(($i+1)) $(($j+1)) $(($k+1)) $(($l+1)) 0 5 8 BFmark"
# else
# 	echo "more than two pairs"
# 	exit 1
# fi
# # prestring="$first $second GREEN Fomark $prestring"
# # echo $prestring
# # exit 0
# # -t 0 means layout mode 0 (default 1)
echo -ne ">$id\n$seq\n$struct" | $VIENNA/bin/RNAplot -t $mode --pre "$prestring " --post "$poststring" # the final space is important to keep "" #"$span GREEN Fomark"

sed -i 's/fsize setlinewidth/9 setlinewidth/' ${id}_ss.ps # change line width

ps2pdf -dEPSCrop ${id}_ss.ps # bounding box

# pdfcrop ${id}_ss.pdf # crop margin automatically
pdfcropmargins -v -u -s ${id}_ss.pdf -o ${id}_ss-crop.pdf # pip install pdfCropMargins
mv ${id}_ss-crop.pdf ${id}.pdf # final output
rm ${id}_ss.p* # remove temp files

echo "output ${id}.pdf"
exit

# code for automatic rotation
# rm -f /tmp/angles
# for angle in `seq 0 30 179` # every 30 degrees
# do
# 	echo "rotating", $angle
# 	pdfjam --angle $angle --outfile /tmp/${id}-${angle}.pdf ${id}.pdf
# 	pdfcrop /tmp/${id}-${angle}.pdf
# 	(echo -ne "$angle\t"; ./bbox.py /tmp/${id}-${angle}-crop.pdf) >> /tmp/angles
# done

# bestangle=`cat /tmp/angles | sort -nk3 | head -1 | cut -f 1` # smallest height (y)
# echo crude best angle, $bestangle

# low=$((bestangle-20))
# high=$((bestangle+20))

# echo trying refined angles from $low to $high

# for angle in `seq $low 5 $high` # every 5 degrees
# do
# 	echo "rotating", $angle
# 	pdfjam --angle $angle --outfile /tmp/${id}-${angle}.pdf ${id}.pdf
# 	pdfcrop /tmp/${id}-${angle}.pdf
# 	(echo -ne "$angle\t"; ./bbox.py /tmp/${id}-${angle}-crop.pdf) >> /tmp/angles
# done

# bestangle=`cat /tmp/angles | sort -nk3 | head -1 | cut -f 1` # smallest height (y)
# echo refined best angle, $bestangle
# rm ${id}.pdf
# mv /tmp/${id}-${bestangle}-crop.pdf ${id}_loop.pdf
# echo "final output: ${id}_loop.pdf"
# echo "------------------------"
# echo