#!/usr/bin/env bash

#FILE:  png2mp4.sh
#DATE:  11 DEC 2023
#AUTH:  G. E. Deschaines
#DESC:  Converts a sequence of PNG files to JPEG files and merges
#       the JPEG files into a MP4 video file.
#
# NOTE:
#
#   [1] Requires ffmpeg and the ImageMagick magick program.
#
#   [2] Invoke this shell script while working within the pyThreeD
#       subdirectory containing img_####.png files (i.e., ./Ximg)
#       as ../util/png2mp4.sh since it's assumed relative path
#       ../pyThreeD.py exists from where this shell script executes.

FFMPEG_EXE=/usr/bin/ffmpeg

# Get list of PNG files.

PNG_LIST=(`ls *.png`)
if [ "${#PNG_LIST[*]}" == "0" ]
then
  echo "error:  No PNG files in present working directory."
  exit -1
fi

# Convert each PNG file into a JPEG file.

JPG_LIST=""
for FILE in ${PNG_LIST[@]}
do
  NAME=${FILE%.png}
  echo "converting:  $FILE"
  convert $FILE $NAME.jpg
  JPG_LIST="$JPG_LIST $NAME.jpg"
done

# Merge all JPEG files into MP4 video file.

FPS=`grep "img_FPS = " ../pyThreeD.py | gawk '{ print gensub(";","",1,$3) }'`

SYSNAM=`uname -a`

if [ "${SYSNAM%% *}" == "Linux" ] && [ -e ${FFMPEG_EXE} ]
then
  echo "Creating MP4 video:  img_anim.mp4 @ 25 fps"
  ${FFMPEG_EXE} -loglevel error -f image2 -framerate ${FPS} -i ./img_%04d.jpg -r 25 -vcodec libx264 -b:v 8000k -crf 18 -pix_fmt yuv420p ./img_anim.mp4
  exit 0
fi

if [ "${SYSNAM%%_*}" == "CYGWIN" ] && [ -e ${FFMPEG_EXE} ]
then
  echo "Creating MP4 video:  img_anim.mp4 @ 25 fps"
  ${FFMPEG_EXE} -loglevel error -f image2 -framerate ${FPS} -i ./img_%04d.jpg -r 25 -vcodec mpeg4 -b:v 8000k -crf 18 -pix_fmt yuv420p ./img_anim.mp4
  exit 0
fi

echo "Could not find ffmpeg to create ./img_anim.mp4 file."
exit 0
