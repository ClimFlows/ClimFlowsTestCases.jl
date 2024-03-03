using ClimFlowsTestCases
using Documenter

DocMeta.setdocmeta!(ClimFlowsTestCases, :DocTestSetup, :(using ClimFlowsTestCases); recursive=true)

makedocs(;
    modules=[ClimFlowsTestCases],
    authors="The ClimFlows contributors",
    sitename="ClimFlowsTestCases.jl",
    format=Documenter.HTML(;
        canonical="https://ClimFlows.github.io/ClimFlowsTestCases.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ClimFLows/ClimFlowsTestCases.jl.git",
    devbranch="main",
)
