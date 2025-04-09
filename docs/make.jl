using Juleanita
using Documenter

DocMeta.setdocmeta!(
        Juleanita, 
        :DocTestSetup, 
        :(using Juleanita); 
        recursive=true)

makedocs(;
    modules=[Juleanita],
    authors="LisaSchlueter <lschlueter@lbl.gov> and contributors",
    sitename="Juleanita.jl",
    format=Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://LisaSchlueter.github.io/Juleanita.jl/stable/",
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
) 

deploydocs(;
    repo="github.com/LisaSchlueter/Juleanita.jl",
    forcepush = true, 
    push_preview = true,
)
