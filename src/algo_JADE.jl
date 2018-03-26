using Distributions

function algo_JADE(f::Function,plotf::Function,vp::VarParam,ap::JADEParam,op::OptParam)
    ϵ_max=1.0e-2
    ϵ_min=1.0e-8
    ϵ=ϵ_max
    maxit=ap.niter
    nind=ap.n
    progr=op.progress
    graph=op.graphics
    log=op.logging
    verb=op.verbose
    v=vp.vars
    nvar=size(v,1)
    p=ap.p
    c=ap.c
    (nobj,neq,nineq)=ftest(f,zeros(nvar))

    μCR=0.5
    μF=0.5
    std=0.1
    CR=ones(Float64,nind)
    F=ones(Float64,nind)
    meanSCR=zeros(Float64)
    meanSF=zeros(Float64)

    if maximum([1,Int(floor(p*nind))])==1
        error("p is such that JADE->DE")
    end

    fopt=zeros(Float64,nobj)
    foptold=zeros(Float64,nobj)
    indopt=Individual(nvar,nobj,neq,nineq)
    indtrueopt=Individual(nvar,nobj,neq,nineq)

    pop=popinitrand(nind,nvar,nobj,neq,nineq,v)

    fpop=SharedArray{Float64}(nobj,nind)
    heqpop=SharedArray{Float64}(neq,nind)
    hineqpop=SharedArray{Float64}(nineq,nind)

    paralleval!(f,pop,fpop,heqpop,hineqpop,log)
    nftot=countfeasible(pop,ϵ)

    maxeq=zeros(neq)
    maxineq=zeros(nineq)
    maxviol!(pop,maxeq,maxineq)

    (imin,itruemin)=popleaderpercent(p,pop,maxeq,maxineq,ϵ)

    copy!(indopt,pop.ind[imin])
    copy!(indtrueopt,pop.ind[itruemin])

    isfeasible(indtrueopt,ϵ) ? fopt[1]=indtrueopt.f[1] : fopt[1]=1.0E38

    popdonor=popinit0(nind,nvar,nobj,neq,nineq)
    poptrial=popinit0(nind,nvar,nobj,neq,nineq)

    plotdata=PlotData()

    #evals=Array{Int64,1}(0)
    #opts=Array{Individual,1}(0)
    #epsvals=Array{Int64,1}(0)
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
        rand!(Cauchy(μF,std),F)
        truncateregen1!(F,μF,std)
        @notime donorgen!(popdonor,pop,indopt,F)

        #perform recombination
        rand!(Normal(μCR,std),CR)
        truncate01!(CR)
        @notime recomb!(poptrial,pop,popdonor,CR)

        #restore box bounds (doesn't alloc)
        @notime fixit!(poptrial,v)

        #evaluate new population
        @notime paralleval!(f,poptrial,fpop,heqpop,hineqpop,log)
        nevals += nind
        @notime nfeasible=countfeasible(poptrial,ϵ)
        nftot += nfeasible

        #find maximum constraint violations (doesn't alloc)
        @notime maxviol!(pop,maxeq,maxineq) #this could be improved by merging with those of previous round

        #greedy selection considering constraints
        @notime greedyselect!(pop,poptrial,maxeq,maxineq,ϵ,CR,F,meanSCR,meanSF)
        if abs(meanSCR[])>0
            μCR=(1.0-c)*μCR+c*meanSCR[]
            μF=(1.0-c)*μF+c*meanSF[]
        end

        #find new best
        @notime begin
        foptold[1]=fopt[1]
        (imin,itruemin)=popleaderpercent(p,pop,maxeq,maxineq,ϵ)
        copy!(indopt,pop.ind[imin])
        copy!(indtrueopt,pop.ind[itruemin])

        isfeasible(indtrueopt,ϵ) ? fopt[1]=indtrueopt.f[1] : fopt[1]=1.0E38
        end

        @notime begin
        if fopt[1]<foptold[1]
            push!(plotdata.evals,nevals)
            logopt=Individual(nvar,nobj,neq,nineq)
            copy!(logopt,indtrueopt)
            push!(plotdata.opts,logopt)
            push!(plotdata.eps,ϵ)
            if progr
                println("It. ",it," New fopt,feasible=",fopt[1]," ",isfeasible(indopt,ϵ))
                #println(indtrueopt)
                if graph plotf(plotdata) end
            end
        end
        end
        ϵ_old=ϵ
        (ϵ,upd)=epsilonupdate(it,maxit,ϵ,ϵ_min,ϵ_max,nfeasible,nind)
        if upd
            #println("eps:",ϵ_old,"->",ϵ)
            fopt[1]=1.0E38
        end
    end #end while
end #end timing

    if graph plotf(plotdata) end
    if verb
        println("JADE on exit μCR,μF,eps=",μCR," ",μF," ",ϵ)
        println("Feasibility:",nftot,"/",nevals,"=",nftot/nevals)
    end
    return indtrueopt
end
