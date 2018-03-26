using Plots
using Parameters

@with_kw struct PlotData
    evals::Array{Int64}=Array{Int64}(0)
    opts::Array{Individual}=Array{Individual}(0)
    eps::Array{Float64}=Array{Float64}(0)
end

function plotinit()
    pyplot(leg=false)
end

function plotend()
    Plots.savefig("final.pdf")
end

function printplot(plotf::Function)
    println("--- Plotting ---")
    println(string(plotf))
end
