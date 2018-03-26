myid()==1 && __precompile__()

module Optopus

  include("util.jl")
  include("variable.jl")
  include("population.jl")
  include("graphics.jl")
  include("algorithm.jl")
  include("optimize.jl")
  include("objective.jl")
  include("constraints.jl")

  include("algo_DE.jl")
  include("algo_JADE.jl")

  export optimize,VParam,OParam,DEParam,JADEParam
  export Variable,Individual,logreader,@nowarn,machinfo,PlotData

#for testing purposes
  export popinitrand,popprint,paralleval!,maxviol!,bubblesort!,bubblesort2!,Population

end
