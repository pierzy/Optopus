function printobj(f::Function)
    println("--- Objective ---")
    println(string(f))
end

function ftest(f::Function,x::Array{Float64})
    (nobj,neq,nineq)=f()
    f1=Vector{Float64}(nobj)
    heq1=Vector{Float64}(neq)
    hineq1=Vector{Float64}(nineq)
    try
        f(x,f1,heq1,hineq1)
    catch
        error("Possible mismatch between objective function contents, please check!")
    end
    return nobj,neq,nineq
end


function paralleval!(f::Function,pop::Population,fpop::SharedArray{Float64,2},heqpop::SharedArray{Float64,2},hineqpop::SharedArray{Float64,2},log::Bool)
    nind=pop.nind[]
    nobj=pop.ind[1].nobj[]
    nvar=pop.ind[1].nvar[]
    neq=pop.ind[1].neq[]
    nineq=pop.ind[1].nineq[]


    #@sync @parallel
    @notime begin
    @inbounds for i=1:nind
        @views f(pop.ind[i].x,fpop[1:nobj,i],heqpop[1:neq,i],hineqpop[1:nineq,i])
    end
    end
    @notime begin
    @inbounds for i=1:nind
        @inbounds for j=1:nobj
            pop.ind[i].f[j]=fpop[j,i]
        end
        @inbounds for j=1:neq
            pop.ind[i].heq[j]=heqpop[j,i]
        end
        @inbounds for j=1:nineq
            pop.ind[i].hineq[j]=hineqpop[j,i]
        end
    end
    end
    if log
        logfile=open("logfile.log","a")
        for i=1:nind
	           writedlm(logfile,[pop.ind[i].x' pop.ind[i].f' pop.ind[i].heq' pop.ind[i].hineq'],"\t")
        end
        close(logfile)
    end
    #error()
    nothing
end

function logreader()
    logfile=open("logfile.log","r")
    data=readdlm(logfile)
    nvars=size(data,2)-1
    close(logfile)
    v=data[:,1:nvars]
    f=data[:,nvars+1:nvars+1]
    return v,f
end
