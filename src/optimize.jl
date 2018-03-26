using StatsBase
using PyPlot
using Parameters

abstract type OptParam end

@with_kw struct OParam <: OptParam
    verbose::Bool=false
    progress::Bool=false
    graphics::Bool=false
    logging::Bool=false
    nrun::Int64=1
    #OParam(verbose)=new(verbose,false,false,false,1)
    #OParam(verbose,progress)=new(verbose,progress,false,false,1)
    #OParam(verbose,progress,graphics)=new(verbose,progress,graphics,false,1)
    #OParam(verbose,progress,graphics,logging)=new(verbose,progress,graphics,logging,1)
    #OParam(verbose,progress,graphics,logging,nrun)=new(verbose,progress,graphics,logging,nrun)
end

Base.show(io::IO, op::OParam) = print(io, "Optimization: verbose=",op.verbose,", progress=",op.progress,", graphics=",op.graphics,", logging=",op.logging,", nrun=",op.nrun)

function optimize(f::Function,plotf::Function,vp::VarParam,ap::AlgParam,op::OptParam)

  verb=op.verbose
  progr=op.progress
  graph=op.graphics
  logging=op.logging

  if verb println("Welcome to Optopus") end

  logging && (logfile=open("logfile.log","w+"); close(logfile))

  checkvars(vp.vars)
  if verb
      printopt(op)
      printvars(vp.vars)
      printobj(f)
      printalg(ap)
  end

  if graph
      verb && printplot(plotf)
      plotinit()
  end

  s=Symbol("algo_$(ap.kind)")
  algo=getfield(Optopus,s)

  runfopt=1.0E38
  allfopt=Array{Float64}(op.nrun)
  runxopt=Array{Float64}(size(vp.vars,1))
  for irun=1:op.nrun
      indopt=algo(f,plotf,vp,ap,op)
      if verb
          println("Done. xopt=",indopt.x," f=",indopt.f[1]," ",indopt.heq," ",indopt.hineq)
          #println("heqopt=",heqopt)
          #println("hineqopt=",hineqopt)
      end
      fopt=indopt.f[1]
      xopt=indopt.x[:]
      allfopt[irun]=fopt
      fopt<runfopt && (runfopt=fopt; runxopt=xopt)
  end

  if verb && op.nrun>1
      println("--- Statistics ---")
      (mean,std)=mean_and_std(allfopt)
      (min,max)=extrema(allfopt)
      println("nrun=",op.nrun," min=",min," avg=",mean," max=",max," std=",std)
  end

  if graph plotend() end

  return runxopt,runfopt

end

function printopt(op::OptParam)
    println("--- Optimization ---")
    println(op)
end
