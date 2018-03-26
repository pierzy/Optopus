using Hwloc

macro nowarn(expr)
    quote
        stderr = STDERR
        stream = open("nul", "a")
        redirect_stderr(stream)
        result = $(esc(expr))
        redirect_stderr(stderr)
        close(stream)
        result
    end
end

macro notime(body)
    quote
        $(esc(body))
    end
end
"""
Views in Julia still allocate some memory (since they need to keep
a reference to the original array). This type allocates no memory
and does no bounds checking. Use it with caution.
"""
mutable struct UnsafeVectorView{Float64} <: AbstractVector{Float64}
    offset::Int64
    len::Int64
    ptr::Ptr{Float64}
end

UnsafeVectorView{Float64}(parent::DenseArray{Float64}, start_ind::Int64, len::Int64) = UnsafeVectorView{Float64}(start_ind - 1, len, pointer(parent))
Base.size(v::UnsafeVectorView) = (v.len,)
Base.getindex(v::UnsafeVectorView, idx) = unsafe_load(v.ptr, idx + v.offset)
Base.setindex!(v::UnsafeVectorView, value, idx) = unsafe_store!(v.ptr, value, idx + v.offset)
Base.length(v::UnsafeVectorView) = v.len
@static if VERSION >= v"0.6-"
    Base.IndexStyle{V <: UnsafeVectorView}(::Type{V}) = Base.IndexLinear()
else
    Base.linearindexing{V <: UnsafeVectorView}(::Type{V}) = Base.LinearFast()
end

"""
UnsafeVectorView only works for isbits types. For other types, we're already
allocating lots of memory elsewhere, so creating a new View is fine.
This function looks type-unstable, but the isbits(T) test can be evaluated
by the compiler, so the result is actually type-stable.
"""
function fastview{Float64}(parent::AbstractArray{Float64}, start_ind::Int64, len::Int64)
        UnsafeVectorView(parent, start_ind, len)
end

function fastreview!(fv::UnsafeVectorView,start_ind::Int64, len::Int64)
    fv.offset=start_ind-1
    fv.len=len
end

"""
Fallback for non-contiguous arrays, for which UnsafeVectorView does not make
sense.
"""
#fastview(parent::AbstractArray, start_ind::Int64, len::Int64) = @view(parent[start_ind:(start_ind + len - 1)])

function machinfo()
topology = Hwloc.topology_load()
summary = Hwloc.getinfo(topology)
println("Machine overview:")
for obj in summary
    obj_type = obj[1]
    count = obj[2]
    println("$count $obj_type")
end
end

function truncate01!(x::Array{Float64})
    for i = eachindex(x)
        if x[i]<0.0
            x[i] = 0.0
        end
        if x[i]>1.0
            x[i] = 1.0
        end
    end
end

function truncateregen1!(x::Array{Float64},μCR::Float64,std::Float64)
    for i = eachindex(x)
        if x[i]<0.0
            x[i]=rand(Cauchy(μCR,std))
        end
        if x[i]>1.0
            x[i] = 1.0
        end
    end
end
