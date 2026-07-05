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
        repo = "https://github.com/biplab37/NJLModels.jl",
        devbranch="main",
        devurl = "dev",
        deploy_url = "https://biplab37.github.io/NJLModels.jl/dev/",
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
