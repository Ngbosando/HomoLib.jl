using Revise
using HomoLib: theorical_bound
using Test
using CairoMakie
using LinearAlgebra


function test_theoretical_bounds_plot()
    f_range = 0.0:0.05:0.9
    E₁ = 1
    E₂ = 50
    ν₁ = 0.45
    ν₂ =0.3
    λ₁ =  E₁*ν₁ / ((1+ν₁)*(1-2ν₁))
    μ₁ = E₁/ 2*(1+ν₁)
    κ₁ = λ₁ + 2*μ₁ / 3
   

    λ₂ =  E₂*ν₂ / ((1+ν₂)*(1-ν₂))
    μ₂ = E₂/ 2*(1+ν₂)
    κ₂ = λ₂ + 2*μ₂ / 3


    k₁, k₂ = 1.0, 5.0

  
        fig = Figure(size = (800, 400))

        ax1 = Axis(fig[1, 1], title = "Elastic bounds κ", xlabel = "Volume Fraction", ylabel = "Bulk Modulus κ")
        ax2 = Axis(fig[1, 2], title = "Thermal bounds k", xlabel = "Volume Fraction", ylabel = "Conductivity k")

        κ_voigt = Float64[]
        κ_reuss = Float64[]
        κ_HS_low = Float64[]
        κ_HS_up = Float64[]
        μ_HS_low = Float64[]
        μ_HS_up = Float64[]

        k_voigt = Float64[]
        k_reuss = Float64[]
        k_HS_low = Float64[]
        k_HS_up = Float64[]

        for f in f_range
            
            push!(κ_voigt, theorical_bound(:elastic, :voigt, κ₁, μ₁, κ₂, μ₂,f).κ)
            push!(κ_reuss, theorical_bound(:elastic, :reuss, κ₁, μ₁, κ₂, μ₂,f).κ)
            hs = theorical_bound(:elastic, :hashin_shtrikman, κ₁, μ₁, κ₂, μ₂,f)
            push!(κ_HS_low, hs.κ_lower)
            push!(κ_HS_up, hs.κ_upper)
        

            push!(k_voigt, theorical_bound(:thermal, :voigt, k₁, k₂, f))
            push!(k_reuss, theorical_bound(:thermal, :reuss, k₁, k₂, f))
            htk = theorical_bound(:thermal, :hashin_shtrikman, k₁, k₂, f)
            push!(k_HS_low, htk.k_lower)
            push!(k_HS_up, htk.k_upper)
        end

        lines!(ax1, f_range, κ_voigt, label = "Voigt")
        lines!(ax1, f_range, κ_reuss, label = "Reuss")
        lines!(ax1, f_range, κ_HS_low, label = "HS Lower")
        lines!(ax1, f_range, κ_HS_up, label = "HS Upper")
        axislegend(ax1)

        lines!(ax2, f_range, k_voigt, label = "Voigt")
        lines!(ax2, f_range, k_reuss, label = "Reuss")
        lines!(ax2, f_range, k_HS_low, label = "HS Lower")
        lines!(ax2, f_range, k_HS_up, label = "HS Upper")
        axislegend(ax2)
        return fig
        
    end



   a = test_theoretical_bounds_plot()

