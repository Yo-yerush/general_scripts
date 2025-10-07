svg_2_tiff <- function(svg_path, dpi = 600, output_file = NULL, width_mm = NULL, width_pixels = NULL) {
    suppressMessages(suppressWarnings(library(rsvg)))
    suppressMessages(suppressWarnings(library(magick)))

    if (is.null(output_file)) {
        if (!grepl("\\.svg$", svg_path, ignore.case = TRUE)) {
            stop("Input file must have .svg extension")
        }
        output_file <- sub("\\.svg$", ".tiff", svg_path, ignore.case = TRUE)
    }


    img <- image_read_svg(svg_path)

    ####################
    # resize by width
    if (!is.null(width_mm) && is.null(width_pixels)) {
        # convert mm to pixels based on dpi
        width_pixels <- width_mm * (dpi / 25.4)
        # width_pixels <- width_mm * 3.7795275591
    }
    if (!is.null(width_pixels)) {
        # resize image to specified width while maintaining aspect ratio
        img <- image_resize(img, paste0(width_pixels, "x"))
    }

    image_write(img, output_file, format = "tiff", density = dpi)
    cat("saved:", basename(output_file), "\n")
}
