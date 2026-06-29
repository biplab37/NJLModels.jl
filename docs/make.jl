using Documenter, NJLModels

makedocs(
    modules  = [NJLModels],
    sitename = "NJLModels",
    authors  = "Biplab Mahato",
    pages    = Any[
        "Home"         => "index.md",
    ]
)

deploydocs(
    repo = "github.com/biplab37/NJLModels.jl.git",
)
