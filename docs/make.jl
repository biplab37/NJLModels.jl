using Documenter, NJLModels, DocumenterVitepress

makedocs(
    modules  = [NJLModels],
    sitename = "NJLModels",
    authors  = "Biplab Mahato",
    doctest=false,
    pages    = [
        "Home"         => "index.md",
    ],
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/biplab37/NJLModels.jl.git",
        devbranch="main",
        devurl = "dev",
    ),
    warnonly = [:missing_docs],
)

DocumenterVitepress.deploydocs(
    repo = "https://github.com/biplab37/NJLModels.jl.git",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
