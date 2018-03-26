abstract type AlgParam end

function printalg(ap::AlgParam)
    println("--- Algorithm ---")
    println(ap)
end

struct DEParam <: AlgParam
    kind::String
    F::Float64
    CR::Float64
    n::Int64
    niter::Int64
    DEParam(F,CR,n,niter)=new("DE",F,CR,n,niter)
end
Base.show(io::IO, ap::DEParam) = print(io, "Algorithm: type=",ap.kind,", F=",ap.F,", CR=",ap.CR,", n=",ap.n,", niter=",ap.niter)

struct JADEParam <: AlgParam
    kind::String
    p::Float64
    c::Float64
    n::Int64
    niter::Int64
    JADEParam(p,c,n,niter)=new("JADE",p,c,n,niter)
end
Base.show(io::IO, ap::JADEParam) = print(io, "Algorithm: type=",ap.kind,", p=",ap.p,", c=",ap.c,", n=",ap.n,", niter=",ap.niter)
