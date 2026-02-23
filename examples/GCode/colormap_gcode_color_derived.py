from openscad import *
from pymachineconfig import *

# read in the machine configuration file.
# Note: if one does not exist, it creates default structures
#   to work from, but does not save the file.  The end user
#   needs to do that explicitly.
mc = MachineConfig()
print("\nConfigfile:",mc.configfile())

# desplay the working parameters that are read in or created
print("\n############")
working = mc.working_config()
print("The default machine set to:")
for k in working.keys():
    print("  ",k,working[k])
print("############")

# Define parameters for the part
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
engrave_color = mc.gen_color2str(power=300,feed=6000)
print("Engrave color:",engrave_color)

text_3d_2 = text_3d
text_3d_2 = text_3d_2.projection(cut=True)
text_3d_2 = text_3d_2.color(engrave_color)
text_3d_2 = text_3d_2.translate([38,50,10])
text_3d_2.show()

# Mark the backplate as cut (power=100% and feed=400mm/min)
cut_color = mc.gen_color2str(power=1000,feed=400)
print("Cut color:",cut_color)

plate_2 = plate
plate_2 = plate_2.projection(cut=True)
plate_2 = plate_2.color(cut_color)
plate_2 = plate_2.translate([0,40,0])
plate_2.show()

