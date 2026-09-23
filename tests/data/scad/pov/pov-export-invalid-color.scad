// Regression for issue #988: unpainted faces (SVG fill:none -> Color4f())
// must fall back to the scheme default in POV export, not emit rgbf <-1,-1,-1,2>.
linear_extrude(height = 5)
    import("../../svg/fill-none-closed.svg", center = true);
