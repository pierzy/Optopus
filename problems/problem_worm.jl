# define variables
nstick=40
v=Array{Variable}(nstick)
low=0.0
high=2*pi
for i=1:nstick
    v[i]=Variable("Real",low,high,"")
end
vp=VParam(v)

#define algorithmic parameters
#ap=DEParam(0.8,0.5,40,2000)
ap=JADEParam(0.05,0.1,40,2000)

#define optimization run parameters
op=OParam(true,true,true,false,1)

#define objective function (2 behaviours)
include("geometry.jl")
function worm!()
    nobj=1
    neq=1
    nineq=0
    return nobj,neq,nineq
end

function worm!(x::Array{Float64},f::AbstractArray{Float64},heq::AbstractArray{Float64},hineq::AbstractArray{Float64})
    n=length(x)
    l=sqrt(2)/n
    p=[0.0,0.0]
    cuts=0.0
    for i=1:n
        o=p
        p += [l*cos(x[i]),l*sin(x[i])]
        l2 = Line(Point{Float64}(o[1], o[2]), Point{Float64}(p[1], p[2]))
        l1 = Line(Point{Float64}(0.4, 0.6), Point{Float64}(0.6, 0.4))
        (inter,dist)=intersection(l1, l2)
        if inter
            cuts += 1.0
            #if dist<=0.2*sqrt(2.0) cuts += dist end
            #if dist>0.2*sqrt(2.0) cuts += dist-0.2*sqrt(2.0) end
        end
    end
    cuts2=p[2]-1.0+p[1]
    dist=sqrt((p[1]-0.51)^2+(p[2]-0.51)^2)
    f[1]=dist
    heq[1]=cuts
    hineq=nothing
    #hineq[2]=-cuts2
    nothing
end

function plotter(x::Vector{Int64},opt::Vector{Individual})
    npt=length(x)
    if npt==0 return end
    y=map(t -> t.f[1], opt)
    ndof=length(opt[1].x)
    l=sqrt(2)/ndof
    px=[0.0]; py=[0.0]
    ppx=0.0; ppy=0.0
    for i=1:ndof
        ppx += l*cos(opt[npt].x[i])
        ppy += l*sin(opt[npt].x[i])
        push!(px,ppx)
        push!(py,ppy)
    end
    p1=plot(x,y,yaxis=:log10)
    p2=plot(px,py,xlims=(0,1),ylims=(0,1))
    p2=scatter!(px,py)
    p2=plot!([0.0,1.0],[0.0,1.0],linecolor=:red)
    p2=plot!([0.4,0.6],[0.6,0.4],linecolor=:black)
    p=plot(p1,p2,layout=(1,2))
    display(p)
end
