svg_2_tiff <- function(svg_path, dpi = 600, output_file = NULL) {
    suppressMessages(suppressWarnings(library(rsvg)))
    suppressMessages(suppressWarnings(library(magick)))

    if (is.null(output_file)) {
        if (!grepl("\\.svg$", svg_path, ignore.case = TRUE)) {
            stop("Input file must have .svg extension")
        }
        output_file <- sub("\\.svg$", ".tiff", svg_path, ignore.case = TRUE)
    }

    img <- image_read_svg(svg_path)
    image_write(img, output_file, format = "tiff", density = dpi)

    cat("saved:", basename(output_file), "\n")
}
