using Documenter
using SpectroscopyTools

makedocs(
    sitename = "SpectroscopyTools.jl",
    modules = [SpectroscopyTools],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Chirp Correction" => "tutorials/chirp_correction.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/garrekstemo/SpectroscopyTools.jl.git",
)
