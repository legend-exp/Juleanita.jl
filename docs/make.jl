using Juleanita
using Documenter

DocMeta.setdocmeta!(Juleanita, :DocTestSetup, :(using Juleanita); recursive=true)

makedocs(;
    modules=[Juleanita],
    authors="LisaSchlueter <lschlueter@lbl.gov> and contributors",
    sitename="Juleanita.jl",
    format=Documenter.HTML(;
        canonical="https://"LisaSchlueter".github.io/Juleanita.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/"LisaSchlueter"/Juleanita.jl",
    devbranch="main",
)
