from openscad import *
import copy

"""
 Colormap methods for mapping object colors to laser power and feed rates. 
"""

# def_Color_Table - the default colormap
_def_Color_Table ={
    "L00": {"power":1.0,"feed":1.0,"color":0x000000},
    "L01": {"power":1.0,"feed":1.0,"color":0x0000FF},
    "L02": {"power":1.0,"feed":1.0,"color":0xFF0000},
    "L03": {"power":1.0,"feed":1.0,"color":0x00E000},
    "L04": {"power":1.0,"feed":1.0,"color":0xD0D000},
    "L05": {"power":1.0,"feed":1.0,"color":0xFF8000},
    "L06": {"power":1.0,"feed":1.0,"color":0x00E0E0},
    "L07": {"power":1.0,"feed":1.0,"color":0xFF00FF},
    "L08": {"power":1.0,"feed":1.0,"color":0xB4B4B4},
    "L09": {"power":1.0,"feed":1.0,"color":0x0000A0},
    "L10": {"power":1.0,"feed":1.0,"color":0xA00000},
    "L11": {"power":1.0,"feed":1.0,"color":0x00A000},
    "L12": {"power":1.0,"feed":1.0,"color":0xA0A000},
    "L13": {"power":1.0,"feed":1.0,"color":0xC08000},
    "L14": {"power":1.0,"feed":1.0,"color":0x00A0FF},
    "L15": {"power":1.0,"feed":1.0,"color":0xA000A0},
    "L16": {"power":1.0,"feed":1.0,"color":0x808080},
    "L17": {"power":1.0,"feed":1.0,"color":0x7D87B9},
    "L18": {"power":1.0,"feed":1.0,"color":0xBB7784},
    "L19": {"power":1.0,"feed":1.0,"color":0x4A6FE3},
    "L20": {"power":1.0,"feed":1.0,"color":0xD33F6A},
    "L21": {"power":1.0,"feed":1.0,"color":0x8CD78C},
    "L22": {"power":1.0,"feed":1.0,"color":0xF0B98D},
    "L23": {"power":1.0,"feed":1.0,"color":0xF6C4E1},
    "L24": {"power":1.0,"feed":1.0,"color":0xFA9ED4},
    "L25": {"power":1.0,"feed":1.0,"color":0x500A78},
    "L26": {"power":1.0,"feed":1.0,"color":0xB45A00},
    "L27": {"power":1.0,"feed":1.0,"color":0x004754},
    "L28": {"power":1.0,"feed":1.0,"color":0x86FA88},
    "L29": {"power":1.0,"feed":1.0,"color":0xFFDB66},
    "T1":  {"power":0.0,"feed":0.0,"color":0xF36926},
    "T2":  {"power":0.0,"feed":0.0,"color":0x0C96D9}
    }

# set the user accessible colormap to the default
Color_Table = copy.deepcopy(_def_Color_Table)

# reset_default_colormap - reset the working colormap to the default
def reset_default_colormap():
    """reset_default_colormap: change potentially modified labeled power,
       feed, and color associations back to their default.  The
       default color table is compatible with LightBurn's 

    """
    global Color_Table
    Color_Table = copy.deepcopy(_def_Color_Table)

# colormap2str - return the working labled color as an OpenSCAD
#   compatible string representation of the hex value starting with a
#   '#'
def colormap2str(key):
    return "#{:X}".format(int(colormap(key)))

# colormap - return the working labled color hex value
def colormap(key):
    return Color_Table[key]["color"]

# powermap - return the working labled power
def powermap(key):
    return Color_Table[key]["power"]

# feedmap - return the working labled feed
def feedmap(key):
    return Color_Table[key]["feed"]

# setpower - overwrite the working labeled power
def set_powermap(key,val):
    Color_Table[key]["power"] = val
    return Color_Table[key]["power"]

# setfeed - overwrite the working labeled feed
def set_feedmap(key,val):
    Color_Table[key]["feed"] = val
    return Color_Table[key]["feed"]

# setcolor - overwrite the working labeled color
def set_colormap(key,val):
    Color_Table[key]["color"] = val
    return Color_Table[key]["color"]

def gen_color(red=-1,green=-1,blue=-1,power=-1,feed=-1):
    if (red!=-1 or green!=-1 or blue!=-1) and (power!=-1 or feed!=-1):
        print("Error (gen_color): can only set either RGB or PF values.")
        raise ValueError("Can only set either RGB or PF values.")
    color = 0
    if red   != -1: color |= (int(255.0*red) << 24)
    if green != -1: color |= (int(255.0*green) << 16)
    if blue  != -1: color |= (int(255.0*blue) << 8)
    if power != -1: color |= (int(255.0*power) << 24)
    if feed  != -1: color |= (int(255.0*feed) << 16)

    return color

def gen_color2str(red=-1,green=-1,blue=-1,power=-1,feed=-1):
    color = gen_color(red,green,blue,power,feed)
    
    return "#{:X}".format(color)
