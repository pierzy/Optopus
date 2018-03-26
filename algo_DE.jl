function algo_DE(f::Function,plotf::Function,vp::VarParam,ap::DEParam,op::OptParam)
    eps=1.0e-6
    maxit=ap.niter
    nind=ap.n
    progr=op.progress
    graph=op.graphics
    log=op.logging
    v=vp.vars
    nvar=size(v,1)
    CR=ones(Float64,nind)*ap.CR
    F=ones(Float64,nind)*ap.F
    meanSCR=zeros(Float64)
    meanSF=zeros(Float64)
    (nobj,neq,nineq)=ftest(f,zeros(nvar))
    fopt=zeros(Float64,nobj)
    foptold=zeros(Float64,nobj)
    indopt=Individual(nvar,nobj,neq,nineq)

    pop=popinitrand(nind,nvar,nobj,neq,nineq,v)

    fpop=SharedArray{Float64}(nobj,nind)
    heqpop=SharedArray{Float64}(neq,nind)
    hineqpop=SharedArray{Float64}(nineq,nind)

    paralleval!(f,pop,fpop,heqpop,hineqpop,log)

    maxeq=zeros(neq)
    maxineq=zeros(nineq)
    maxviol!(pop,maxeq,maxineq)

    imin=popleader(pop,maxeq,maxineq,eps)

    copy!(indopt,pop.ind[imin])

    isfeasible(indopt,eps) ? fopt[1]=pop.ind[imin].f[1] : fopt[1]=1.0E38

    popdonor=popinit0(nind,nvar,nobj,neq,nineq)
    poptrial=popinit0(nind,nvar,nobj,neq,nineq)

    evals=Array{Int64,1}(0)
    opts=Array{Individual,1}(0)
    nevals=nind

    foptold[1]=1.0E38
    #apFm=-ap.F
    it=0
    winner=zeros(Int64)
    wfeasible=zeros(Bool)

@notime begin
    while it<maxit
        it += 1
        #generate donors by mutation
        @notime donorgen!(popdonor,pop,indopt,F)

        #perform recombination
        @notime recomb!(poptrial,pop,popdonor,CR)

        #restore box bounds (doesn't alloc)
        @notime fixit!(poptrial,v)

        #evaluate new population
        @notime paralleval!(f,poptrial,fpop,heqpop,hineqpop,log)
        nevals += nind

        #find maximum constraint violations (doesn't alloc)
        @notime maxviol!(pop,maxeq,maxineq) #this could be improved by merging with those of previous round

        #greedy selection considering constraints
        @notime greedyselect!(pop,poptrial,maxeq,maxineq,eps,CR,F,meanSCR,meanSF)

        #find new best
        @notime begin
        foptold[1]=fopt[1]
        imin=popleader(pop,maxeq,maxineq,eps)
        copy!(indopt,pop.ind[imin])
        isfeasible(indopt,eps) ? fopt[1]=indopt.f[1] : fopt[1]=1.0E38
        end

        @notime begin
        if fopt[1]<foptold[1]
            push!(evals,nevals)
            logopt=Individual(nvar,nobj,neq,nineq)
            copy!(logopt,indopt)
            push!(opts,logopt)
            if progr
                println("It. ",it," New fopt,feasible=",fopt[1]," ",isfeasible(indopt,eps))
                #println(indopt)
                if graph plotf(evals,opts) end
            end
        end

        end
    end
end
    if graph
        plotf(evals,opts)
    end

    return indopt
end

function recomb!(poptrial::Population,pop::Population,popdonor::Population,CR::Vector{Float64})
nind=poptrial.nind[]
nvar=poptrial.ind[1].nvar[]
iir=zeros(Int64,1)
rr=zeros(Float64,nvar)
for j=1:nind
    rand!(iir,1:nvar)
    rand!(rr)
    for k=1:nvar
        #(rr[k]<ap.CR || ir[1]==k) ? poptrial.ind[j].x[k]=popdonor.ind[j].x[k] : poptrial.ind[j].x[k]=pop.ind[j].x[k]
        if rr[k]<CR[j]
            poptrial.ind[j].x[k]=popdonor.ind[j].x[k]
        elseif iir[1]==k
            poptrial.ind[j].x[k]=popdonor.ind[j].x[k]
        else
         poptrial.ind[j].x[k]=pop.ind[j].x[k]
        end
    end
end
end

function donorgen!(popdonor::Population, pop::Population,indopt::Individual,F::Vector{Float64})
nind=popdonor.nind[]
nvar=popdonor.ind[1].nvar[]
ir=zeros(Int64,3)
for j=1:nind
    rand3!(ir,nind) #other option is usin StatsBase ir=sample(1:nind, 3, replace=false)
#        istrat=rand(1:2)
#        if istrat==1

        #DE/best/1/bin
        BLAS.copy!(popdonor.ind[j].x,1,indopt.x,1,nvar)
        BLAS.axpy!(nvar,F[j],pop.ind[ir[2]].x,1,popdonor.ind[j].x,1)
        BLAS.axpy!(nvar,-F[j],pop.ind[ir[3]].x,1,popdonor.ind[j].x,1)

        #DE/rand/1/bin
        #BLAS.copy!(popdonor.ind[j].x,1,pop.ind[ir[1]].x,1,nvar)
        #BLAS.axpy!(nvar,F[j],pop.ind[ir[2]].x,1,popdonor.ind[j].x,1)
        #BLAS.axpy!(nvar,-F[j],pop.ind[ir[3]].x,1,popdonor.ind[j].x,1)

        #for k in 1:nvar
            #popdonor.ind[j].x[k]=indopt.x[k]
        #    popdonor.ind[j].x[k] += ap.F.*pop.ind[ir[2]].x[k]
        #    popdonor.ind[j].x[k] += -ap.F.*pop.ind[ir[3]].x[k]
        #end
#        else
#            popdonor.ind[j].x[1:nvar]=pop.ind[ir[1]].x[1:nvar]+ap.F*(pop.ind[ir[2]].x[1:nvar]-pop.ind[ir[3]].x[1:nvar])
#        end
end
end

function greedyselect!(pop::Population,poptrial::Population,maxeq::Vector{Float64},maxineq::Vector{Float64},eps::Float64,CR::Vector{Float64},F::Vector{Float64},meanSCR::Array{Float64},meanSF::Array{Float64})
const nind=pop.nind[]
meanSCR[]=0.0
meanSF[]=0.0
nsucc::Int64=0
num::Float64=0.0
den::Float64=0.0
for j=1:nind
    winner=epsiloncompare(pop.ind[j],poptrial.ind[j],maxeq,maxineq,eps)
    if winner==2
        copy!(pop.ind[j],poptrial.ind[j])
        nsucc += 1
        meanSCR[] += CR[j]
        num += F[j]^2
        den += F[j]
    end
    #winner==2 && swapprint(pop.ind[j],poptrial.ind[j])
end
if nsucc>0
     meanSCR[]=meanSCR[]/nsucc
     meanSF[]=num/den
end
end

function rand3!(ir::Vector{Int64},irange::Int64)
    rand!(ir,1:irange)
    while ir[1]==ir[2] || ir[1]==ir[2] || ir[2]==ir[3]
        rand!(ir,1:irange)
    end
end
