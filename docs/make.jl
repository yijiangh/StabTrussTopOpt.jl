using Documenter, StabTrussTopOpt

makedocs(;
    modules=[StabTrussTopOpt],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/yijiangh/StabTrussTopOpt.jl/blob/{commit}{path}#L{line}",
    sitename="StabTrussTopOpt.jl",
    authors="Yijiang Huang",
    assets=[],
)

deploydocs(;
    repo="github.com/yijiangh/StabTrussTopOpt.jl",
)
