import Base:copy!

struct Individual
    nvar::Array{Int64,0}
    nobj::Array{Int64,0}
    neq::Array{Int64,0}
    nineq::Array{Int64,0}
    x::Array{Float64,1}
    f::Array{Float64,1}
    heq::Array{Float64,1}
    hineq::Array{Float64,1}
    viol::Array{Float64,0}
    Individual(nv::Int64,no::Int64,ne::Int64,ni::Int64)=new(nv*ones(),no*ones(),ne*ones(),ni*ones(),zeros(nv),zeros(no),zeros(ne),zeros(ni),zeros())
    Individual(nv::Int64,no::Int64,ne::Int64,ni::Int64,x::Array{Float64,1},f::Array{Float64,1},heq::Array{Float64,1},hineq::Array{Float64,1},viol::Float64)=new(nv*ones(),no*ones(),ne*ones(),ni*ones(),x,f,heq,hineq,viol)
end
Base.copy!(dest::Individual,src::Individual)=  for k=1:9 copy!(getfield(dest,k),getfield(src,k)) end #change 9 to something else if Individual changes

function isfeasible(a::Individual,eps::Float64)
    feasible=true
        any(x->abs(x)-eps>0.0,a.heq) && (feasible=false)

        any(x->x>0.0,a.hineq) && (feasible=false)
    return feasible
end

struct Population
    nind::Int64
    ind::Array{Individual}
    Population(nind)=new(nind,Array{Individual,1}(nind))
end

function countfeasible(pop::Population,eps::Float64)
    nf=0
    nind=nind=pop.nind[]
    for i=1:nind
        isfeasible(pop.ind[i],eps) && (nf += 1)
    end
    return nf
end

function popinitrand(nind::Int64,nvar::Int64,nobj::Int64,neq::Int64,nineq::Int,v::Array{Variable})
    thepop=Population(nind)
    rr=zeros(Float64,nvar)
    for ipop=1:nind
        thepop.ind[ipop]=Individual(nvar,nobj,neq,nineq)
        rand!(rr)
        for ivar=1:nvar
            thepop.ind[ipop].x[ivar]=v[ivar].min+rr[ivar]*(v[ivar].max-v[ivar].min)
        end
    end
    return thepop
end

function popinit0(nind::Int64,nvar::Int64,nobj::Int64,neq::Int64,nineq::Int64)
    thepop=Population(nind)
    for ipop=1:nind
        thepop.ind[ipop]=Individual(nvar,nobj,neq,nineq)
        for ivar=1:nvar
            thepop.ind[ipop].x[ivar]=0.0
        end
    end
    return thepop
end

function maxviol!(pop::Population,maxeq::Array{Float64},maxineq::Array{Float64})
    nind=pop.nind[]
    neq=pop.ind[1].neq[]
    nineq=pop.ind[1].nineq[]

    maxeq[:]=-1.0E38
    maxineq[:]=-1.0E38
    for i=1:nind
        for j=1:neq
            abs(pop.ind[i].heq[j])>maxeq[j] && (maxeq[j]=abs(pop.ind[i].heq[j]))
        end
        for j=1:nineq
            pop.ind[i].hineq[j]>maxineq[j] && (maxineq[j]=pop.ind[i].hineq[j])
        end
    end

    for j=1:neq
        maxeq[j]<0.0 && (maxeq[j]=1.0)
        #maxeq[j]=1.0
    end
    for j=1:nineq
        maxineq[j]<0.0 && (maxineq[j]=1.0)
        #maxineq[j]=1.0
    end

    nothing
end


function popleader(pop::Population,maxeq::Array{Float64,1},maxineq::Array{Float64,1},eps::Float64)

    nind=pop.nind[]
    neq=pop.ind[1].neq[]
    nineq=pop.ind[1].nineq[]

    imin=1

    for i=2:nind
        winner=epsiloncompare(pop.ind[imin],pop.ind[i],maxeq,maxineq,eps)
        winner==2 && (imin=i)
    end

    return imin
end

function popleaderpercent(p::Float64,pop::Population,maxeq::Array{Float64,1},maxineq::Array{Float64,1},eps::Float64)
    nind=pop.nind[]
    iorder=collect(1:nind)
    bubblesort2!(iorder,pop,maxeq,maxineq,eps)
    top=rand(1:maximum([1,Int(floor(p*nind))]))
    imin=iorder[top]
    itruemin=iorder[1]
    return imin,itruemin
end

function epsiloncompare(a::Individual,b::Individual,maxeq::Array{Float64},maxineq::Array{Float64},eps::Float64)
    fa=isfeasible(a,eps)
    fb=isfeasible(b,eps)

    (fa && fb)  &&  (a.f[1]<b.f[1] ? winner=1 : winner=2;  return winner) #both are feasible, return best

    (fa && !fb) && (winner=1;  return winner)

    (fb && !fa) && (winner=2;  return winner)

    consa=0.0
    consb=0.0
    nobj=a.nobj[]
    neq=a.neq[]
    nineq=a.nineq[]

    if neq>0
        consa=consa+fsumvscal!(a.heq,maxeq,neq)
        consb=consb+fsumvscal!(b.heq,maxeq,neq)
    end

    if nineq>0
        consa=consa+fsumvscal!(a.hineq,maxineq,nineq)
        consb=consb+fsumvscal!(b.hineq,maxineq,nineq)
    end

    consa<consb ? winner=1 : winner=2
    return winner

end

function fsumvscal!(c::Array{Float64,1},d::Array{Float64,1},ne::Int64)
    s=0.0
    for i in 1:ne
        s += abs(c[i])/d[i]
    end
    return s
end

function popprint(pop::Population,descr::String)
    println(descr)
    for j=1:pop.nind
        println(pop.ind[j].x," ",pop.ind[j].f[1]," ",pop.ind[j].heq)
    end
end

function swapprint(a::Individual,b::Individual)
    println(a.x," ",a.f," ",a.heq," ",a.hineq,"->",b.x," ",b.f," ",b.heq," ",b.hineq,)
end

function bubblesort!(pop::Population,maxeq::Array{Float64},maxineq::Array{Float64},eps::Float64)
    passes::Int64 = 0
    nvar=pop.ind[1].nvar[]
    nobj=pop.ind[1].nobj[]
    neq=pop.ind[1].neq[]
    nineq=pop.ind[1].nineq[]
    tmp=Individual(nvar,nobj,neq,nineq)
    while(true)
        no_swaps = bubbleswap!(tmp,pop, passes,maxeq,maxineq,eps)
        if no_swaps == 0
            break
        end
        passes += 1
    end
    return pop
end

function bubbleswap!(tmp::Individual,pop::Population, passes::Int64, maxeq::Vector{Float64},maxineq::Vector{Float64},eps::Float64)
    no_swaps::Int64 = 0
    nind=pop.nind[]
    for vec_index in 1:(nind - 1 - passes)
        winner=epsiloncompare(pop.ind[vec_index],pop.ind[vec_index + 1],maxeq,maxineq,eps)
        if winner==2
            no_swaps += 1
            copy!(tmp,pop.ind[vec_index]) #dest source
            copy!(pop.ind[vec_index],pop.ind[vec_index + 1])
            copy!(pop.ind[vec_index + 1],tmp)
        end
    end
    return no_swaps
end

function bubblesort2!(iorder::Vector{Int64},pop::Population,maxeq::Array{Float64},maxineq::Array{Float64},eps::Float64)
    passes::Int64 = 0
    while(true)
        no_swaps = bubbleswap2!(iorder,pop,passes,maxeq,maxineq,eps)
        if no_swaps == 0
            break
        end
        passes += 1
    end
    return iorder
end

function bubbleswap2!(iorder::Vector{Int64},pop::Population, passes::Int64, maxeq::Vector{Float64},maxineq::Vector{Float64},eps::Float64)
    no_swaps::Int64 = 0
    nind=pop.nind[]
    for vec_index in 1:(nind - 1 - passes)
        winner=epsiloncompare(pop.ind[iorder[vec_index]],pop.ind[iorder[vec_index + 1]],maxeq,maxineq,eps)
        if winner==2
            no_swaps += 1
            tmp=iorder[vec_index]
            iorder[vec_index]=iorder[vec_index + 1]
            iorder[vec_index + 1]=tmp
        end
    end
    return no_swaps
end
