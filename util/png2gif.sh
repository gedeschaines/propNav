#!/usr/bin/env bash

#FILE:  png2gif.sh
#DATE:  11 DEC 2023
#AUTH:  G. E. Deschaines
#DESC:  Converts a sequence of PNG files to GIF files and merges
#       the GIF files into an animated GIF file.
#
# NOTE:
#
#   [1] Requires ImageMagick magick program.
#
#   [2] Invoke this shell script while working within the pyThreeD
#       subdirectory containing img_####.png files (i.e., ./Ximg)
#       as ../util/png2gif.sh since it's assumed relative path
#       ../pyThreeD.py exists from where this shell script executes.

# Get list of PNG files.

PNG_LIST=(`ls *.png`)
if [ "${#PNG_LIST[*]}" == "0" ]
then
  echo "error:  No PNG files in present working directory."
  exit -1
fi

# Convert each PNG file into a GIF file.

GIF_LIST=""
for FILE in ${PNG_LIST[@]}
do
  NAME=${FILE%.png}
  echo "converting:  $FILE"
  convert $FILE $NAME.gif
  GIF_LIST="$GIF_LIST $NAME.gif"
done

# Calculate image delay time in hundredths of a second.

FPS=`grep "img_FPS = " ../pyThreeD.py | gawk '{ print gensub(";","",1,$3) }'`
DELAY=`echo "100 $FPS" | gawk '{ print $1/$2 }'`

# Merge all GIF files into the animated gif file.

# NOTE: The -size WxH pixel values in the convert expression below must
#       match that specified by fig = plt.figure(figsize=(w,h), dpi=###)
#       statement in ../threeD.py script, where W = w*dpi and H = h*dpi.

echo "Creating animated gif file:  img_anim.gif @ 0.0$DELAY sec/image"
convert -size 800x600 -dispose None -delay $DELAY $GIF_LIST -loop 2 img_anim.gif

exit 0

