# define variables
v=Array{Variable}(2)
low=-1.2
high=1.1
v[1]=Variable("Real",low,high,"x1")
v[2]=Variable("Real",low,high,"x2")
vp=VParam(v)

#define algorithmic parameters
ap=JADEParam(0.05,0.1,40,8000)
#ap=DEParam(0.8,0.5,40,8000)

#define optimization run parameters
op=OParam(true,false,false,false,50)

#define objective function (2 behaviours)
function g11!()
    nobj=1
    neq=1
    nineq=0
    return nobj,neq,nineq
end

function g11!(x::Array{Float64},f::AbstractArray{Float64},heq::AbstractArray{Float64},hineq::AbstractArray{Float64})
    heq[1]=x[2]-x[1]^2
    f[1]=x[1]^2+(x[2]-1.0)^2
    hineq=nothing
    nothing
end
