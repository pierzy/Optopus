push!(LOAD_PATH, pwd())
#addprocs(4)
@everywhere using Optopus
@everywhere include("problem_g03.jl") #load problem definition
#machinfo()

#srand(1234)
#(xopt,fopt)=optimize(g11!,vp,ap,op)
#for --track-allocation
#Profile.clear_malloc_data()
@time (xopt,fopt)=optimize(g03!,plotter,vp,ap,op)

#using Traceur
#@trace (xopt,fopt)=optimize(g08!,vp,ap,op)
#(allv,allf)=logreader()
