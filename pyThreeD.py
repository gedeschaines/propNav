# -*- coding: utf-8 -*-
"""
FILE:  pyThreeD.py
DATE:  06 DEC 2023
AUTH:  G. E. Deschaines
DESC:  Three dimensional (3D) drawing of objects defined as collections
       of polygons.
REFS:
     
  [1] This Python script was refactored from threeD.c available at:
      https://github.com/gedeschaines/threeD/blob/master/src/threeD.c


Disclaimer:  See DISCLAIMER file.

"""

import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("* Error: Matplotlib package required.")
    print("         Suggest installing the SciPy stack.")
    sys.exit()


global fig       # Instantiated matplotlib pyplot figure
global draw3D    # Instantiated Draw3D object
global started   # Draw3D.MainLoop method started flag

from draw3D import Draw3D

started = False


def onClick(event):
    """
    Mouse button pressed handler
    """
    global fig
    global draw3D
    global started
    
    #b = event.button
    #x = event.xdata
    #y = event.ydata
    if not started:
       started = True
       draw3D.MainLoop()
       if draw3D.doneflag:
           fig.clear()
           plt.close(fig)
           sys.exit()
       started = False
       

if __name__ == "__main__":
    
    global fig
    global draw3D
    
    # Set image save flag and rate (FPS).
    
    img_Save = False  # Save rendered images as PNG files
    img_FPS = 50;     # (semi-colon is hold over from threeD.c)
    
    # Process command line arguments.
    
    CaseId = None  # TXYZ.OUT.CaseId
    
    if len(sys.argv) > 2:
        CaseId   = sys.argv[1]
        img_Save = False if sys.argv[2] == "0" else True
    else:
        print("usage:  python3 pyThreeD.py CaseId [0|1]")
        print("where:  CaseId - read trajectory data from ./out/TXYZ.OUT.CaseId")
        print("        [0|1]  - save rendered image flag (0=False, 1=True)")
        sys.exit()
        
    File  = "./out/TXYZ.OUT." + CaseId
    title = "pyThreeD - " + File
    
    # Instantiate a matplotlib pyplot figure.
    
    #fig = plt.figure(figsize=(6.0,6.0), dpi=100.0)     # 465x462 viewport
    fig = plt.figure(figsize=(8.0,6.0), dpi=100.0)     # 620x462 viewport
    #fig = plt.figure(figsize=(8.0,8.0), dpi=100.0)     # 620x616 viewport
    #fig = plt.figure(figsize=(7.74,7.78), dpi=100.0)   # 600x599 viewport
    #fig = plt.figure(figsize=(10.0,8.0), dpi=100.0)    # 775x616 viewport
    #fig = plt.figure(figsize=(10.32,7.80), dpi=100.0)  # 800x601 viewport
    
    # Specify plotting layout and parameters.
    
    ax = fig.add_subplot(111, autoscale_on=False, animated=False)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    bbox = ax.bbox.get_points()
    w = round(bbox[1,0] - bbox[0,0])
    h = round(bbox[1,1] - bbox[0,1])
    #print("w= %d, h= %d" % (w, h))
    ax.set_xlim(0.0, w)
    ax.set_ylim(h, 0.0)
    ax.set_facecolor('#d8dcd6')  # light gray
    fig.canvas.draw()
    
    # Instantiate a Draw3D object.
    
    draw3D = Draw3D(txyzFile=File, w=w, h=h, fig=fig, ax=ax,
                    imgSave=img_Save, imgFPS=img_FPS)
    
    # Load object shape polygon data.
    
    draw3D.ReadPolyData()
    
    # Assign mouse button press handler.
    
    cidbtn = fig.canvas.mpl_connect('button_press_event', onClick)
    
    # Assign key press event handler.
    
    cidkey = fig.canvas.mpl_connect('key_press_event', draw3D.onPress)
    
    # Present instructions, then display figure window for rendering. 
    
    print("With cursor in the displayed Figure, press a mouse button")
    print("to initiate threeD animation. If animation is exited or has")
    print("completed, close figure to terminate this program or press")
    print("a mouse button to restart animation. During threeD program")
    print("animation the following key presses are recognized:\n")
    
    print("Press t key to toggle field-of-view towards target.")
    print("Press m key to toggle field-of-view towards missile.")
    print("Press h key to toggle field-of-view along missile heading.")
    print("Press z key to reset zoom to one.")
    print("Press up arrow key to increase zoom.")
    print("Press down arrow key to decrease zoom.")
    print("Press left arrow key to slow animation down by 50 msec increments.")
    print("Press right arrow key to speed animation up by 50 msec increments.")
    print("Press Space key to toggle pause/unpause.")
    print("Press x key to exit animation (press a mouse button to restart).")
    print("Press Esc key to exit program.")
    
    plt.show(block=True)
    
    # Execution terminated, delete draw3D object and exit.
    
    del draw3D
 