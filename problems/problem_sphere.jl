# define variables
v=Array{Variable}(2)
low=-10.0
high=10.0
v[1]=Variable("Real",low,high,"x1")
v[2]=Variable("Real",low,high,"x2")
vp=VParam(v)

#define algorithmic parameters
ap=DEParam(0.8,0.5,20,40)

#define optimization run parameters
op=OParam(true,false,true,false,1)

#define objective function (2 behaviours)
function sphere!()
    nobj=1
    neq=0
    nineq=0
    return nobj,neq,nineq
end

function sphere!(x::Array{Float64},f::AbstractArray{Float64},heq::AbstractArray{Float64},hineq::AbstractArray{Float64})
    heq=nothing
    f[1]=(x[1]-1.0)^2+(x[2]-1.0)^2
    hineq=nothing
    #sleep(1)
    nothing
end

function plotter(x::Vector{Int64},opt::Vector{Individual})
    n=length(x)
    y=Array{Float64}(n)
    for i=1:n
        y[i]=opt[i].f[]
    end
    theplot=plot(x,y,yaxis=:log10)
    display(theplot)
end
