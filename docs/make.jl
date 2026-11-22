using NonmetricMDS
using Documenter

DocMeta.setdocmeta!(NonmetricMDS, :DocTestSetup, :(using NonmetricMDS); recursive=true)

makedocs(;
    modules=[NonmetricMDS],
    checkdocs=:exports,
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/NonmetricMDS.jl/blob/{commit}{path}#{line}",
    sitename="NonmetricMDS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/NonmetricMDS.jl",
        repolink="https://github.com/EvoArt/NonmetricMDS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/NonmetricMDS.jl",
)
