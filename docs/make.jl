using Literate
using Documenter

# Process the literate file first
LITERATE_DIR = joinpath(@__DIR__, "src")
GENERATED_DIR = joinpath(@__DIR__, "src/generated")
isdir(GENERATED_DIR) && rm(GENERATED_DIR, recursive=true)
mkpath(GENERATED_DIR)

Literate.markdown(
    joinpath(LITERATE_DIR, "sensitivity_analysis.jl"),
    GENERATED_DIR,
    documenter=true,
)

# Then build docs with minimal configuration
makedocs(
    sitename = "FEM Analysis of Composites",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Numerical Analysis" => "generated/sensitivity_analysis.md",
    ]
)