@echo off

rem FILE:  png2mp4.bat
rem DATE:  17 DEC 2023
rem AUTH:  G. E. Deschaines
rem DESC:  Converts a sequence of PNG files to JPEG files and merges
rem        the JPEG files into a MP4 video file.
rem
rem NOTE:
rem
rem   [1] Requires ffmpeg and the ImageMagick magick program.
rem
rem   [2] Invoke this batch file while working within the pyThreeD
rem       subdirectory containing img_####.png files (i.e., ./Ximg)
rem       as ../util/png2mp4.bat since it's assumed relative path
rem       ../pyThreeD.py exists from where this batch file executes.

SetLocal EnableExtensions EnableDelayedExpansion

set FFMPEG_EXE=C:\bin\ffmpeg.exe

rem Get list of PNG files.

set JPG_LIST=
set njpg=0
set count=0
for /F "tokens=*" %%f in ( 'dir /B *.png' ) do (
  set /a count=!count! + 1
  set fpng=%%f
  set PNG_LIST[!count!]=!fpng!
  set fjpg=!fpng:~0,-4!.jpg
  if exist !fjpg! (
    set JPG_LIST=!JPG_LIST! !fjpg!
    set /a njpg=!njpg! + 1
  )
)
echo "!njpg!"

rem Check if PNG to JPEG conversion can be skipped.
rem NOTE: Assumes each JPEG file was previously converted
rem       from like named PNG file using this script.

if !count! GTR 0 (
  if !njpg! EQU !count! goto find_fps
)

if !count! EQU 0 (
  rem Get list of JPEG files.
  set JPG_LIST=
  set njpg=0
  for /F "tokens=*" %%f in ( 'dir /B *.jpg' ) do (
    set /a njpg=!njpg! + 1
    set JPG_LIST=!JPG_LIST! %%f
  )
  if !njpg! EQU 0 (
     echo "error:  No PNG or JPEG files in present working directory."
     goto exit0
  )
  rem PNG to JPEG conversion can be skipped.
  goto find_fps
)

rem Convert each PNG file into a JPEG file.

set JPG_LIST=
for /L %%i in (1, 1, !count!) do (
  set NAME=!PNG_LIST[%%i]:~0,-4!
  echo "converting:  !PNG_LIST[%%i]!"
  magick !PNG_LIST[%%i]! !NAME!.jpg
  set JPG_LIST=!JPG_LIST! !NAME!.jpg
)

rem Merge all JPEG files into MP4 video file.

:find_fps
set n=0
for /F "tokens=*" %%t in ( 'find "img_FPS = " ..\pyThreeD.py' ) do (
  set /a n=!n! + 1
  set line=%%t
  if !n! EQU 2 (
     set FPS=!line:~10,2!
     goto create_mp4
  )
)
set FPS=50

:create_mp4
if exist !FFMPEG_EXE! (
  echo "Creating MP4 video:  img_anim.mp4 @ 25 fps"
  !FFMPEG_EXE! -loglevel error -f image2 -framerate !FPS! -i .\img_%%04d.jpg -r 25 -vcodec libx264 -b:v 8000k -crf 18 -pix_fmt yuv420p .\img_anim.mp4
  goto exit0
)

echo "Could not find ffmpeg to create ./img_anim.mp4 file."

:exit0

EndLocal

exit /b 0
