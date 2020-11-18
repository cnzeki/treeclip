#include "svg.h"


const std::string SvgBase::svg_xml_start[] =
{ "<?xml version=\"1.0\" standalone=\"no\"?>\n"
"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"\n"
"\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\n"
"<svg width=\"",
"\" height=\"",
"\" viewBox=\"0 0 ",
"\" version=\"1.0\" xmlns=\"http://www.w3.org/2000/svg\">\n\n"
};
const std::string SvgBase::path_end_poly[] =
{ "\"\n style=\"fill:",
"; fill-opacity:",
"; fill-rule:",
"; stroke:",
"; stroke-opacity:",
"; stroke-width:",
";\"/>\n\n"
};
const std::string SvgBase::path_end_line[] =
{ "\"\n style=\"fill:none; stroke:",
"; stroke-opacity:",
"; stroke-width:",
";\"/>\n\n"
};