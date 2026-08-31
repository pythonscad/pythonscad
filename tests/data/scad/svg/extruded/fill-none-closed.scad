// Regression: closed SVG paths with fill:none must extrude as opaque solids.
// Previously the fill:none paint was stored as alpha=0, so side walls vanished
// and only boolean-cap faces remained (a flat sheet that flipped with the camera).
difference() {
    linear_extrude(height = 15)
        import("../../../svg/fill-none-closed.svg", center = true);
    translate([-30, -20, -10]) cube([60, 40, 10.001]);
    translate([-30, -20, 14.999]) cube([60, 40, 10]);
}
