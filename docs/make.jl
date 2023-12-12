using Documenter, AstroRepresentations
using Pkg 

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    authors="Andrea Pasquale",
    sitename="AstroRepresentations.jl",
    modules=[AstroRepresentations],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Trasformations" => "transform.md"
    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/AstroRepresentations.jl", branch="gh-pages"
)