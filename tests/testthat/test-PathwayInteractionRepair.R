test_that("network labels retain source semantics and undirected edges have no arrow", {
    topology <- list(protEdges = data.frame(
        src = "a", dest = c("b", "c"), directed = c(1L, 2L), reaction_type = c(4L, 2L),
        source_reaction_type = c("Control(Out: INHIBITION of BiochemicalReaction)", "Binding")
    ), mixedEdges = data.frame(src = "c", dest = "d", reaction_type = 6L))
    constants <- new.env(parent = emptyenv())
    utils::data("reaction_type", package = "SpaMTP", envir = constants)
    reaction <- constants$reaction_type
    edges <- SpaMTP:::.pn_prepare_edges(topology, reaction)
    expect_identical(edges$reaction_name[edges$dest == "b"], topology$protEdges$source_reaction_type[1])
    expect_identical(edges$arrowhead[edges$dest == "b"], "thead")
    expect_identical(edges$arrowhead[edges$dest == "c"], "none")
    expect_identical(edges$reaction_name[edges$dest == "d"], "Process(activation)")
})
