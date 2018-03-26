push!(LOAD_PATH, pwd())
#addprocs(4)
@everywhere using Optopus
@everywhere include("problem_g11.jl") #load problem definition
#machinfo()

epsil=1.0E-6
nind=5
nobj=1
neq=1
nineq=0
pop=popinitrand(nind,2,nobj,neq,nineq,vp.vars)
fpop=SharedArray{Float64}(nobj,nind)
heqpop=SharedArray{Float64}(neq,nind)
hineqpop=SharedArray{Float64}(nineq,nind)

paralleval!(g11!,pop,fpop,heqpop,hineqpop)
popprint(pop,"")
maxeq=zeros(neq)
maxineq=zeros(nineq)
maxviol!(pop,maxeq,maxineq)
iorder=collect(1:nind)
@time bubblesort2!(iorder,pop,maxeq,maxineq,epsil)
popprint(pop,"")

srand(1234)
#(xopt,fopt)=optimize(g11!,vp,ap,op)
#for --track-allocation
#Profile.clear_malloc_data()
#@time (xopt,fopt)=optimize(g11!,vp,ap,op)

#using Traceur
#@trace (xopt,fopt)=optimize(g08!,vp,ap,op)
#(allv,allf)=logreader()
