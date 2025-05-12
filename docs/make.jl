using Juleanita
using Documenter
using Literate 

DocMeta.setdocmeta!(
        Juleanita, 
        :DocTestSetup, 
        :(using Juleanita); 
        recursive=true)

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end

gen_content_dir = joinpath(@__DIR__, "src", "tutorials")
for tut_lit_fn in filter(fn -> endswith(fn, "_lit.jl"), readdir(gen_content_dir))
    lit_src_fn = joinpath(gen_content_dir, tut_lit_fn)
    tut_basename = tut_lit_fn[1:end-7] # remove "_lit.jl"
    Literate.notebook(lit_src_fn, gen_content_dir, name = tut_basename, documenter = true, credit = true, execute = false)
    Literate.markdown(lit_src_fn, gen_content_dir, name = tut_basename, documenter = true, credit = true, postprocess = fix_literate_output)
end

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
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Basic I/O" => "tutorials/reading_data.md",
        ],
        "API" => "api.md",
         "ML Quality Cuts" => "ML_QC.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = false,
    warnonly = ("nonstrict" in ARGS),
) 

deploydocs(;
    repo="github.com/LisaSchlueter/Juleanita.jl",
    forcepush = true, 
    push_preview = true,
)
