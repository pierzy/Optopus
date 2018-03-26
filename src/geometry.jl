struct Point{T}
    x::T
    y::T
end

struct Line{T}
    s::Point{T}
    e::Point{T}
end

function intersection(l1::Line{T}, l2::Line{T}) where T<:Real
    a1 = l1.e.y - l1.s.y
    b1 = l1.s.x - l1.e.x
    c1 = a1 * l1.s.x + b1 * l1.s.y

    a2 = l2.e.y - l2.s.y
    b2 = l2.s.x - l2.e.x
    c2 = a2 * l2.s.x + b2 * l2.s.y

    Δ = a1 * b2 - a2 * b1
    # If lines are parallel, intersection point will contain infinite values
    inter=Point((b2 * c1 - b1 * c2) / Δ, (a1 * c2 - a2 * c1) / Δ)
    if isnan(inter.x) return false,0.0 end
    if isnan(inter.y) return false,0.0 end
    if isinf(inter.x) return false,0.0 end
    if isinf(inter.y) return false,0.0 end

    v1=[l1.e.x-l1.s.x,l1.e.y-l1.s.y]
    v2=[inter.x-l1.s.x,inter.y-l1.s.y]
    proj=dot(v2,v1/norm(v1))
    if proj>norm(v1) || proj<0.0 return false,0.0 end
    return true,proj
end
