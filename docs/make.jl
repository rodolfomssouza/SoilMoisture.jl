using SoilMoisture
using Documenter

DocMeta.setdocmeta!(SoilMoisture, :DocTestSetup, :(using SoilMoisture); recursive=true)

makedocs(;
    modules=[SoilMoisture],
    authors="Rodolfo Souza rodolfomssouza@gmail.com",
    repo="https://github.com/rodolfomssouza/SoilMoisture.jl/blob/{commit}{path}#{line}",
    sitename="SoilMoisture.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://rodolfomssouza.github.io/SoilMoisture.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/rodolfomssouza/SoilMoisture.jl",
    devbranch="main",
)
