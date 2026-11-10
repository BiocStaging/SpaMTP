# Development assistance: see inst/COPYRIGHTS.

#' Embed a pathway network in an HTML document
#'
#' Includes the interactive viewer created by [pathwayNetworkPlots()] in an
#' HTML report or vignette. The viewer and its data are copied into the
#' document, so the original file is no longer needed after rendering.
#' A download button saves the standalone network HTML from the embedded copy.
#'
#' @param file HTML file returned by [pathwayNetworkPlots()].
#' @param height Preferred viewer height in pixels; a positive finite number.
#'   The viewer adapts to the browser height, with a minimum of 420 pixels
#'   or the requested height when smaller.
#' @param title Accessible title for the embedded viewer.
#' @return An `htmltools` tag containing the viewer and its download button.
#'   Print it in an HTML R Markdown chunk to display the network. Other output
#'   formats should show the `network_summary` attribute of `file` instead.
#' @details The embedded viewer uses its bundled JavaScript and runs offline.
#'   Its styles are isolated from the surrounding document in an iframe. Only
#'   embed HTML from a trusted source; JavaScript in the supplied file executes
#'   when the document is viewed. This function does not launch a browser or
#'   require a browser during document builds.
#' @seealso [pathwayNetworkPlots()]
#' @examples
#' # Small offline HTML for illustrating the embedding interface.
#' file <- tempfile(fileext = ".html")
#' writeLines("<!doctype html><html><body>Network preview</body></html>", file)
#' viewer <- embedPathwayNetwork(file, height = 400)
#' unlink(file)
#' # The complete HTML is retained even after the input file is removed.
#' inherits(viewer, "shiny.tag")
#' @export
embedPathwayNetwork <- function(
    file, height = 720,
    title = "Interactive pathway network"
) {
    if (!is.character(file) || length(file) != 1L || is.na(file) ||
        !file.exists(file) || dir.exists(file)) {
        stop("file must name an existing pathway network HTML file.",
            call. = FALSE
        )
    }
    if (!is.numeric(height) || length(height) != 1L || !is.finite(height) ||
        height <= 0) {
        stop("height must be one positive finite number.", call. = FALSE)
    }
    if (!is.character(title) || length(title) != 1L || is.na(title) ||
        !nzchar(trimws(title))) {
        stop("title must be one non-empty string.", call. = FALSE)
    }
    document <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"),
        collapse = "\n"
    )
    if (!nzchar(trimws(document))) {
        stop("The network HTML file is empty.", call. = FALSE)
    }
    download <- paste(
        "var frame=this.closest('.spamtp-network-embed').querySelector('iframe');",
        "var url=window.URL.createObjectURL(new window.Blob([frame.srcdoc],",
        "{type:'text/html;charset=utf-8'}));",
        "var link=document.createElement('a');",
        "link.href=url;link.download='spamtp-pathway-network.html';",
        "document.body.appendChild(link);link.click();link.remove();",
        "setTimeout(function(){window.URL.revokeObjectURL(url);},1000);"
    )
    htmltools::tags$div(
        class = "spamtp-network-embed",
        htmltools::tags$p(htmltools::tags$button(
            type = "button", class = "btn btn-default", onclick = download,
            "Download network HTML"
        )),
        htmltools::tags$iframe(
            srcdoc = document, title = title, loading = "lazy",
            style = paste0(
                "width:100%;height:", format(height, scientific = FALSE),
                "px;max-height:90vh;min-height:min(",
                format(height, scientific = FALSE),
                "px,420px);border:1px solid #d7dde5;border-radius:10px;"
            )
        )
    )
}
