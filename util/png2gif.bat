@echo off

rem FILE:  png2gif.bat
rem DATE:  17 DEC 2023
rem AUTH:  G. E. Deschaines
rem DESC:  Converts a sequence of PNG files to GIF files and merges
rem        the GIF files into an animated GIF file.

rem NOTE:  Requires ImageMagick magick program.

SetLocal EnableExtensions EnableDelayedExpansion

rem Get list of PNG files.

set GIF_LIST=
set ngif=0
set count=0
for /F "tokens=*" %%f in ( 'dir /B *.png' ) do (
  set /a count=!count! + 1
  set fpng=%%f
  set PNG_LIST[!count!]=!fpng!
  set fgif=!fpng:~0,-4!.gif
  if exist !fgif! (
     set GIF_LIST=!GIF_LIST! !fgif!
     set /a ngif=!ngif! + 1
  )
)

rem Check if PNG to GIF conversion can be skipped.
rem NOTE: Assumes each GIF file was previously converted
rem       from like named PNG file using this script.
if !count! GTR 0 (
  if !ngif! EQU !count! goto find_fps
)

if !count! EQU 0 (
  rem Get list of GIF files.
  set GIF_LIST=
  set ngif=0
  for /F "tokens=*" %%f in ( 'dir /B *.gif' ) do (
    set /a ngif=!ngif! + 1
    set GIF_LIST=!GIF_LIST! %%f
  )
  if !ngif! EQU 0 (
     echo "error:  No PNG or GIF files in present working directory."
     goto exit0
  )
  rem PNG to GIF conversion can be skipped.
  goto find_fps
)

rem Convert each PNG file into a GIF file.

set GIF_LIST=
for /L %%i in (1, 1, !count!) do (
  set NAME=!PNG_LIST[%%i]:~0,-4!
  echo "converting:  !PNG_LIST[%%i]!"
  magick !PNG_LIST[%%i]! !NAME!.gif
  set GIF_LIST=!GIF_LIST! !NAME!.gif
)

rem Calculate image delay time in hundredths of a second.

:find_fps
set n=0
for /F "tokens=*" %%t in ( 'find "img_FPS = " ..\pyThreeD.py' ) do (
  set /a n=!n! + 1
  set line=%%t
  if !n! EQU 2 (
     set FPS=!line:~10,2!
     set /a DELAY=100 / !FPS!
     goto create_anim_gif
  )
)
set FPS=50
set DELAY=2

rem Merge all GIF files into the animated gif file.
rem NOTE: The -size WxH pixel values in the magick expression below must
rem       match that specified by fig = plt.figure(figsize=(w,h), dpi=###)

:create_anim_gif
echo "Creating animated gif file:  img_anim.gif @ 0.0!DELAY! sec/image"
magick -size 800x600 -dispose None -delay !DELAY! !GIF_LIST! -loop 2 img_anim.gif

:exit0

EndLocal

exit /b 0
