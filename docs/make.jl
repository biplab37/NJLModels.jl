using Documenter, NJLModels, DocumenterVitepress

makedocs(
    modules=[NJLModels],
    sitename="NJLModels.jl",
    authors="Biplab Mahato",
    doctest=false,
    repo=Remotes.GitHub("biplab37", "NJLModels.jl"),
    pages=[
        "Home" => "index.md",
        "Parameters" => "pages/parameters.md",
        "Massgap" => "pages/massgap.md",
        "Indices" => "pages/indices.md",
    ],
    format=DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/biplab37/NJLModels.jl",
        devbranch="main",
        # devurl="dev",
        # deploy_url="https://biplab37.github.io/NJLModels.jl/dev/",
    ),
    warnonly=[:missing_docs, :doctest],
)

DocumenterVitepress.deploydocs(
    repo="github.com/biplab37/NJLModels.jl",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
