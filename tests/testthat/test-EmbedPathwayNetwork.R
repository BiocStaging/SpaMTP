# Development assistance: see inst/COPYRIGHTS.

test_that("an embedded network remains usable after its temporary input is removed", {
    file <- tempfile(fileext = ".html")
    document <- paste0('<!doctype html><html><body><h1>A &amp; B</h1>',
        '<script>window.networkLabel="A & B";</script></body></html>')
    writeLines(document, file)
    viewer <- embedPathwayNetwork(file, height = 480, title = 'A "quoted" network')
    unlink(file)

    frame <- viewer$children[[2L]]
    expect_identical(frame$name, "iframe")
    expect_identical(frame$attribs$srcdoc, document)
    expect_identical(frame$attribs$title, 'A "quoted" network')
    rendered <- htmltools::renderTags(viewer)$html
    expect_match(rendered, "&lt;!doctype html&gt;", fixed = TRUE)
    expect_match(rendered, "A &quot;quoted&quot; network", fixed = TRUE)
    expect_false(grepl(file, rendered, fixed = TRUE))
    expect_match(rendered, "Download network HTML", fixed = TRUE)
    expect_match(rendered, "height:480px", fixed = TRUE)
})

test_that("network embedding rejects unusable inputs before rendering", {
    expect_error(embedPathwayNetwork(tempfile()), "existing pathway")
    expect_error(embedPathwayNetwork(tempdir()), "existing pathway")
    file <- tempfile(fileext = ".html")
    withr::defer(unlink(file))
    writeLines("", file)
    expect_error(embedPathwayNetwork(file), "empty")
    writeLines("<!doctype html><html><body>Network</body></html>", file)
    expect_error(embedPathwayNetwork(file, height = NA_real_), "positive finite")
    expect_error(embedPathwayNetwork(file, height = Inf), "positive finite")
    expect_error(embedPathwayNetwork(file, height = 0), "positive finite")
    expect_error(embedPathwayNetwork(file, title = ""), "non-empty")
})
