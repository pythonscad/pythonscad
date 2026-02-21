from openscad import *
from pycolormap import *

# Define parameters
base_width = 75
base_height = 20
base_thickness = 5
engrave_depth = 1
text_string = "PythonSCAD"
font_size = 8

# Create the Base Plate
plate = cube([base_width, base_height, base_thickness])

# Create the Text (using linear_extrude to give it 
text_3d = text(text_string, size=font_size, font="Arial:style=Bold", halign="center", valign="center").linear_extrude(height=base_thickness)

# Position the text
# We need to move the text into the center of the plate and slightly
# above the bottom surface to avoid rendering artifacts (z-fighting)
positioned_text = text_3d.translate([base_width/2, base_height/2, base_thickness - engrave_depth])

# Perform the Engraving (Difference)
#   Subtract the text shape from the base plate
engraved_object = plate - positioned_text

# Render in 3D (commented out to test ExportGCode
#engraved_object.show()

# Since the GCode generated for the laser cutter does not utilize
# depth (z), offset the engraved portions so that it is not masked by
# the cut parts, but it will be rendered flat on the plane.

#############################
# Use the direct power/feed method.

# Mark the text as engraved (power=50% and feed=100% max)
engrave_color = gen_color2str(power=500,feed=20000)
print("Engrave at:",engrave_color)

text_3d_2 = text_3d
text_3d_2 = text_3d_2.projection(cut=True)
text_3d_2 = text_3d_2.color(engrave_color)
text_3d_2 = text_3d_2.translate([38,50,10])
text_3d_2.show()

# Mark the backplate as cut (power=100% and feed=400mm/min)
cut_color = gen_color2str(power=1000,feed=400)
print("Cut at:",cut_color)

plate_2 = plate
plate_2 = plate_2.projection(cut=True)
plate_2 = plate_2.color(cut_color)
plate_2 = plate_2.translate([0,40,0])
plate_2.show()

#############################
# Use the color-table power/feed method.
#
# modify the color-table entries for cut and engrave
#  cut:
set_powermap("L01",900)
set_feedmap("L01",4000)
#  engrave
set_powermap("L02",350)
set_feedmap("L02",6000)

# engrave the second part
engrave_color = gen_color2str(power=powermap("L02"),feed=feedmap("L02"))
text_3d_3 = text_3d
text_3d_3 = text_3d_3.projection(cut=True)
text_3d_3 = text_3d_3.color(engrave_color)
text_3d_3 = text_3d_3.translate([38,80,10])
text_3d_3.show()

# cut the second part
cut_color = gen_color2str(power=powermap("L01"),feed=feedmap("L01"))
plate_3 = plate
plate_3 = plate_3.projection(cut=True)
plate_3 = plate_3.color(cut_color)
plate_3 = plate_3.translate([0,70,0])
plate_3.show()

print("Second Engrave at:",engrave_color)
print("Second Cut at:",cut_color)
