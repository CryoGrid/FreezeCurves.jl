using FreezeCurves
using Documenter

IS_LOCAL = haskey(ENV,"LOCALDOCS") && ENV["LOCALDOCS"] == "true"

const modules = [
    FreezeCurves,
]

makedocs(modules=modules,
         sitename="FreezeCurves.jl",
         authors="Brian Groenke, Moritz Langer, Jan Nitzbon",
         format=Documenter.HTML(prettyurls=!IS_LOCAL),
         pages=["Home" => "index.md",
                "Inference" => "inference.md",
                "API reference" => [
                       "FreezeCurves" => "api/FreezeCurves.md",
                ],
                "Contributing" => "contributing.md",
])

deploydocs(repo="github.com/CryoGrid/FreezeCurves.jl.git", push_preview=true)
