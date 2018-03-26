# define variables
v=Array{Variable}(2)
low=0.0
high=10.0
v[1]=Variable("Real",low,high,"x1")
v[2]=Variable("Real",low,high,"x2")
vp=VParam(v)

#define algorithmic parameters
#ap=DEParam(0.8,0.5,40,5000)
ap=JADEParam(0.05,0.1,40,5000)

#define optimization run parameters
op=OParam(true,false,false,false,50)


#define objective function (2 behaviours)
function g08!()
    nobj=1
    neq=0
    nineq=2
    return nobj,neq,nineq
end

function g08!(x::Array{Float64},f::AbstractArray{Float64},heq::AbstractArray{Float64},hineq::AbstractArray{Float64})
     f[1]=-(sin(2.0*pi*x[1])^3*sin(2*pi*x[2]))/(x[1]^3*(x[1]+x[2]))
     heq=nothing
     hineq[1]=x[1]^2-x[2]+1.0
     hineq[2]=1.0-x[1]+(x[2]-4.0)^2
     nothing
end
