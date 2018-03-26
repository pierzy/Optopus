# define variables
v=Array{Variable}(10)
n=10
low=0.0
high=1.0
for i=1:n
    v[i]=Variable("Real",low,high,"")
end
vp=VParam(v)

#define algorithmic parameters
ap=JADEParam(0.05,0.1,40,8000)
#ap=DEParam(0.8,0.5,40,8000)

#define optimization run parameters
op=OParam(verbose=true,progress=false,graphics=true,logging=false,nrun=30)

#define objective function (2 behaviours)
function g03!()
    nobj=1
    neq=1
    nineq=0
    return nobj,neq,nineq
end

function g03!(x::Array{Float64},f::AbstractArray{Float64},heq::AbstractArray{Float64},hineq::AbstractArray{Float64})
    n=10
    nd=10.0
    prd=1.0
    s=0.0
    for i=1:n
        prd *= x[i]
        s += x[i]^2
    end
    heq[1]=s-1.0
    f[1]=-sqrt(nd)^n*prd
    hineq=nothing
    nothing
end

function plotter(plotdata::PlotData)
    x=plotdata.evals
    opt=plotdata.opts
    eps=plotdata.eps
    npt=length(x)
    if npt==0 return end
    y1=map(t -> t.f[1]+1.2, opt)  #for log plotting
    println("Points plotted:",npt)
    p1=plot(x,y1,yaxis=:log10)
    p1=plot!(twinx(),x,abs.(eps),color=:red,yaxis=:log10)
    display(p1)
end
