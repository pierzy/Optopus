function fixit!(pop::Population,v::Array{Variable})
    nind=pop.nind
    nvar=pop.ind[1].nvar[]
    for i=1:nind
        for j=1:nvar
            pop.ind[i].x[j]<v[j].min && (pop.ind[i].x[j]=(v[j].min+v[j].max)/2.0)
            pop.ind[i].x[j]>v[j].max && (pop.ind[i].x[j]=(v[j].min+v[j].max)/2.0)
        end
    end
end

function epsilonupdate(it::Int64,maxit::Int64,eps::Float64,epsmin::Float64,epsmax::Float64,nfeasible::Int64,nind::Int64)
    epsold=eps
    progress=it/maxit
    kickin=0.7
    #goldenratio=1.61803398875
    goldenratio=1.5
    #if progress>kickin
    #    eps=epsmin
    #else
    #    eps=(kickin*maxit-it)/(kickin*maxit)*epsmax+epsmin
    #end
    if nfeasible<0.05*nind
        eps=eps*goldenratio
    elseif nfeasible>0.5*nind
        eps=eps/goldenratio
    end
    eps<epsmin && (eps=epsmin)
    eps>epsmax && (eps=epsmax)
    progress>kickin && (eps=epsmin)

    if abs(epsold-eps)<1.0e-9
        upd=false
    else
        upd=true
    end
    return eps,upd
end
