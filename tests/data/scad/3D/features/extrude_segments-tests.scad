extrude()
{
    square(10);
    rotate([0,0,90])
    translate([0,0,50])
    square(10);
}

translate([20,0,0])
extrude(segments=8)
{
    square(10);
    rotate([0,0,90])
    translate([0,0,50])
    square(10);
}

translate([40,0,0])
extrude(segments=32)
{
    square(10);
    rotate([0,0,90])
    translate([0,0,50])
    square(10);
}