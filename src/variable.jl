abstract type VarParam end

struct Variable
  kind::String
  min::Float64
  max::Float64
  descr::String
  Variable(kind,min,max,descr)=new(kind,min,max,descr)
  Variable(kind,min,max)=new(kind,min,max,"")
end

Base.show(io::IO, v::Variable) = print(io, "Var: type=",v.kind,", min=",v.min,", max=",v.max,", descr=",v.descr)

struct VParam <: VarParam
    vars::Array{Variable}
    VParam(vars)=new(vars)
end

function checkvars(v::Array{Variable})
  nvars=size(v,1)
  allowed=Dict("Real"=>1)
  for i=1:nvars
    if !(haskey(allowed,v[i].kind))
      error("Variable type not allowed:",v[i].kind)
    end
  end
end

function printvars(v::Array{Variable})
    nvars=size(v,1)
    println("--- Variables ---")
    for i=1:nvars
        println(i," ",v[i])
    end
end
