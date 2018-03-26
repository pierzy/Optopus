# define variables
v=Array{Variable}(2)
v[1]=Variable("Real",13.0,100.0,"x1")
v[2]=Variable("Real",0.0,100.0,"x2")
vp=VParam(v)

#define algorithmic parameters
ap=JADEParam(0.05,0.1,40,8000)
#ap=DEParam(0.8,0.5,40,8000)

#define optimization run parameters
op=OParam(true,false,false,false,50)

#define objective function (2 behaviours)
function g06!()
    nobj=1
    neq=0
    nineq=2
    return nobj,neq,nineq
end

function g06!(x::Array{Float64},f::AbstractArray{Float64},heq::AbstractArray{Float64},hineq::AbstractArray{Float64})
    hineq[1]=-(x[1]-5.0)^2-(x[2]-5)^2+100.0
    hineq[2]=(x[1]-6.0)^2+(x[2]-5)^2-82.81
    f[1]=(x[1]-10.0)^3+(x[2]-20.0)^3
    heq=nothing
    nothing
end
