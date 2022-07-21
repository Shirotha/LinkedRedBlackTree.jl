include("header.jl")
using .LinkedRedBlackTree

using Luxor

struct Vec2{T <: Number} <: Number
    x::T
    y::T
end
Base.:+(a::Vec2{T}, b::Vec2{T}) where T = Vec2{T}(a.x + b.x, a.y + b.y)
Base.:-(a::Vec2{T}, b::Vec2{T}) where T  = Vec2{T}(a.x - b.x, a.y - b.y)
Base.:*(a::Vec2{T}, b::Number) where T = Vec2{T}(a.x * b, a.y * b)
Base.:*(a::Number, b::Vec2) = b * a
Base.:/(a::Vec2{T}, b::Number) where T = Vec2{T}(a.x / b, a.y / b)
ortho(p::Vec2{T}) where T = Vec2{T}(-p.y, p.x)
det(a::Vec2{T}, b::Vec2{T}) where T = a.x * b.y - a.y * b.x
norm(p::Vec2) = sqrt(p.x^2 + p.y^2)
distance(a::Vec2{T}, b::Vec2{T}) where T = norm(b - a)
Base.length(p::Vec2) = 2
Base.iterate(p::Vec2) = (p.x, 1)
Base.iterate(p::Vec2, state) = state >= 2 ? nothing : (p.y, 2)
Base.min(a::Vec2{T}, b::Vec2{T}) where T = Vec2(min(a.x, b.x), min(a.y, b.y))
Base.max(a::Vec2{T}, b::Vec2{T}) where T = Vec2(max(a.x, b.x), max(a.y, b.y))

Vec2{T}(p::Vec2{S}) where {T, S} = Vec2{T}(p.x, p.y)
Base.convert(::Type{Vec2{T}}, p::Vec2) where T = Vec2{T}(p)
Base.promote_rule(::Type{Vec2{T}}, ::Type{Vec2{S}}) where {T, S} = Vec2{promote_type(T, S)}


abstract type AbstractFace{T} end

mutable struct Site{T, F <: AbstractFace{T}}
    index::Int
    point::Vec2{T}
    weight::T
    face::Union{Nothing, F}

    Site(index, p::Vec2{T}, weight::T, face) where T = new{T, Face{T}}(index, p, weight, face)
end

mutable struct Vertex{T}
    point::Vec2{T}
end

mutable struct HalfEdge{T, F <: AbstractFace{T}}
    origin::Union{Nothing, Vertex{T}}
    destination::Union{Nothing, Vertex{T}}
    twin::Union{Nothing, HalfEdge{T, F}}
    parent::F
    prev::Union{Nothing, HalfEdge{T, F}}
    next::Union{Nothing, HalfEdge{T, F}}

    HalfEdge(o, d, twin, p::AbstractFace{T}, prev, next) where T = new{T, Face{T}}(o, d, twin, p, prev, next)
end

mutable struct Face{T} <: AbstractFace{T}
    site::Site{T}
    children::Union{Nothing, HalfEdge{T, Face{T}}}
end

mutable struct Voronoi{T}
    sites::Vector{Site{T}}
    faces::Vector{Face{T}}
    vertices::Vector{Vertex{T}}
    halfedges::Vector{HalfEdge{T, Face{T}}}

    function Voronoi(points::AbstractVector{Tuple{Vec2{T}, T}}) where T
        sites = Vector{Site{T}}(undef, length(points))
        faces = Vector{Face{T}}(undef, length(points))
        @inbounds for (i, p) in enumerate(points)
            sites[i] = Site(i, p..., nothing)
            faces[i] = Face(sites[i], nothing)
            sites[i].face = faces[i]
        end

        new{T}(sites, faces, Vector{Vertex{T}}(), Vector{HalfEdge{T, Face{T}}}())
    end
end
function Base.push!(v::Voronoi{T}, p::Vec2{T}) where T
    ver = Vertex(p)
    push!(v.vertices, ver)
    return ver
end
function Base.push!(v::Voronoi{T}, f::Face{T}) where T
    he = HalfEdge(nothing, nothing, nothing, f, nothing, nothing)
    push!(v.halfedges, he)
    isnothing(f.children) && (f.children = he)
    return he
end
function Base.delete!(v::Voronoi{T}, ver::Vertex{T}) where T
    i = findfirst(x -> x === ver, v.vertices)
    !isnothing(i) && deleteat!(v.vertices, i)
end
function Base.delete!(v::Voronoi{T}, he::HalfEdge{T}) where T
    i = findfirst(x -> x === he, v.halfedges)
    !isnothing(i) && deleteat!(v.halfedges, i)
end

abstract type AbstractEvent{T} end

mutable struct Arc{T, E <: AbstractEvent{T}}
    site::Site{T}
    left::Union{Nothing, HalfEdge{T, Face{T}}}
    right::Union{Nothing, HalfEdge{T, Face{T}}}
    event::Union{Nothing, E}

    Arc(site::Site{T}, l, r, e) where T = new{T, Event{T}}(site, l, r, e)
end
Arc(site) = Arc(site, nothing, nothing, nothing)

mutable struct Event{T} <: AbstractEvent{T}
    type::Symbol
    y::T
    index::Int

    # type == :site
    site::Union{Nothing, Site{T}}

    Event(site::Site{T}) where T = new{T}(:site, site.point.y, -1, site, nothing, nothing)
    # type == :circle
    point::Union{Nothing, Vec2{T}}
    arc::Union{Nothing, Arc{T}}

    Event(y::T, point::Vec2{T}, arc::Arc{T}) where T = new{T}(:circle, y, -1, nothing, point, arc)
end

mutable struct AABB{T}
    min::Vec2{T}
    max::Vec2{T}
    AABB(a::Vec2{T}, b::Vec2{T}) where T = new{T}(min(a, b), max(a, b))
end
AABB(a::Vec2{T}, b::Vec2{S}) where {T, S} = AABB(promote(a, b)...)
AABB(bounds::AABB) = AABB(bounds.min, bounds.max)
AABB(xmin::T, ymin::T, xmax::T, ymax::T) where T = AABB(Vec2(xmin, ymin), Vec2(xmax, ymax))
AABB(xmin, ymin, xmax, ymax) = AABB(promote(xmin, ymin, xmax, ymax)...)
AABB{T}(xmin, ymin, xmax, ymax) where T = AABB(convert.(T, (xmin, ymin, xmax, ymax))...)
AABB{T}() where T = AABB(zeros(T, 4)...)

topleft(bounds::AABB) = Vec2(bounds.min.x, bounds.max.y)
bottomleft(bounds::AABB) = Vec2(bounds.min.x, bounds.min.y)
bottomright(bounds::AABB) = Vec2(bounds.max.x, bounds.min.y)
topright(bounds::AABB) = Vec2(bounds.max.x, bounds.max.y)

function Base.push!(v::Voronoi{T}, bounds::AABB{T}, side::Symbol) where T
    if side == :left
        return push!(v, topleft(bounds))
    elseif side == :bottom
        return push!(v, bottomleft(bounds))
    elseif side == :right
        return push!(v, bottomright(bounds))
    elseif side == :top
        return push!(v, topright(bounds))
    else
        error("unknown side $side")
    end
end

function expand!(bounds::AABB{T}, point::Vec2{T}) where T
    bounds.min = min(bounds.min, point)
    bounds.max = max(bounds.max, point)
end
isinside(bounds::AABB{T}, point::Vec2{T}) where T =
    bounds.min.x - eps(T) <= point.x <= bounds.max.x + eps(T) &&
    bounds.min.y - eps(T) <= point.y <= bounds.max.y + eps(T)

SIDES = [:left, :bottom, :right, :top]
mutable struct RaycastHit{T}
    point::Vec2{T}
    side::Symbol
end
RaycastHit{T}() where T = RaycastHit(Vec2{T}(0, 0), :none)
function raycast(bounds::AABB{T}, origin::Vec2{T}, direction::Vec2{T}) where T
    t::T = Inf
    result = RaycastHit{T}()
    if direction.x > zero(T)
        t = (bounds.max.x - origin.x) / direction.x
        result.side = :right
        result.point = origin + t * direction
    elseif direction.x < zero(T)
        t = (bounds.min.x - origin.x) / direction.x
        result.side = :left
        result.point = origin + t * direction
    end

    if direction.y > zero(T)
        newT = (bounds.max.y - origin.y) / direction.y
        if newT < t
            result.side = :top
            result.point = origin + newT * direction
        end
    elseif direction.y < zero(T)
        newT = (bounds.min.y - origin.y) / direction.y
        if newT < t
            result.side = :bottom
            result.point = origin + newT * direction
        end
    end
    return result
end
function raycastall(bounds::AABB{T}, origin::Vec2{T}, destination::Vec2{T}) where T
    direction = destination - origin
    result = Vector{Tuple{T, RaycastHit{T}}}()

    if origin.x < bounds.min.x - eps(T) || destination.x < bounds.min.x - eps(T)
        t = (bounds.min.x - origin.x) / direction.x
        if eps(T) < t < one(T) - eps(T)
            p = origin + t * direction
            if bounds.min.y - eps(T) <= p.y <= bounds.max.y + eps(T)
                hit = RaycastHit(p, :left)
                push!(result, (t, hit))
            end
        end
    end

    if origin.x > bounds.max.x + eps(T) || destination.x > bounds.max.x + eps(T)
        t = (bounds.max.x - origin.x) / direction.x
        if eps(T) < t < one(T) - eps(T)
            p = origin + t * direction
            if bounds.min.y - eps(T) <= p.y <= bounds.max.y + eps(T)
                hit = RaycastHit(p, :right)
                push!(result, (t, hit))
            end
        end
    end

    if origin.y < bounds.min.y - eps(T) || destination.y < bounds.min.y - eps(T)
        t = (bounds.min.y - origin.y) / direction.y
        if eps(T) < t < one(T) - eps(T)
            p = origin + t * direction
            if bounds.min.x - eps(T) <= p.x <= bounds.max.x + eps(T)
                hit = RaycastHit(p, :bottom)
                push!(result, (t, hit))
            end
        end
    end

    if origin.y > bounds.max.y + eps(T) || destination.y > bounds.max.y + eps(T)
        t = (bounds.max.y - origin.y) / direction.y
        if eps(T) < t < one(T) - eps(T)
            p = origin + t * direction
            if bounds.min.x - eps(T) <= p.x <= bounds.max.x + eps(T)
                hit = RaycastHit(p, :top)
                push!(result, (t, hit))
            end
        end
    end

    return map(last, sort(result, by = first))
end

mutable struct LinkedVertex{T}
    vertex::Vertex{T}
    prev::Union{Nothing, HalfEdge{T, Face{T}}}
    next::Union{Nothing, HalfEdge{T, Face{T}}}
end

fortune(points::AbstractVector{Vec2{T}}, bounds::AABB{T} = AABB{T}()) where T = fortune(collect(zip(points, zeros(T, length(points)))), bounds)
# TODO: make weights actually do something
function fortune(points::AbstractVector{Tuple{Vec2{T}, T}}, bounds::AABB{T} = AABB{T}()) where T <: Real
    bounds = AABB(bounds)
    voronoi = Voronoi(points)

    events = LRBTree{Event{T}}(; by = e -> e.y, rev = true)
    append!(events, Event.(voronoi.sites))

    y = zero(T)
    root = (a, b) -> begin
        d1 = 1 / 2(a.y - y)
        d2 = 1 / 2(b.y - y)
        @assert !isinf(d1) && !isinf(d2)
        α = d1 - d2
        β = 2(d2 * b.x - d1 * a.x)
        γ = (a.x^2 + a.y^2 - y^2)d1 - (b.x^2 + b.y^2 - y^2)d2
        Δ = β^2 - 4 * α * γ
        return (sqrt(Δ) - β) / 2α
    end

    add = (a_, b_, c_) -> begin
        a, b, c = a_.site.point, b_.site.point, c_.site.point
        # computeConvergencePoint(a, b, c)
        v1 = ortho(a - b)
        v2 = ortho(b - c)
        Δ = (c - a) / 2
        t = det(Δ, v2) / det(v1, v2)
        center = (a + b) / 2 + t * v1
        r = distance(center, a)
        _y = center.y - r

        # isMovingRight(a, b)
        leftmovingright = a.y < b.y
        # isMovingRight(b, c)
        rightmovingright = b.y < c.y
        # getInitialX(a, b, leftmovingright)
        leftx = leftmovingright ? a.x : b.x
        # getInitialX(b, c, rightmovingright)
        rightx = rightmovingright ? b.x : c.x

        isbelow = _y <= y
        isvalid = ((leftmovingright && leftx < center.x) ||
                  (!leftmovingright && leftx > center.x)) &&
                  ((rightmovingright && rightx < center.x) ||
                  (!rightmovingright && rightx > center.x))

        if isbelow && isvalid
            e = Event(_y, center, b_)
            b_.event = e
            insert!(events, e)
        end
    end
    beachline = LRBTree{Arc{T, Event{T}}}()
    rel = @inline (x) -> let d = x
        (
            @inline((a, b) -> d < root(a.site.point, b.site.point)),
            @inline((a, b) -> d > root(a.site.point, b.site.point))
        )
    end

    while !isempty(events)
        event = pop!(events)
        y = event.y
        if event.type == :site
            if isempty(beachline)
                insert!(beachline, Arc(event.site))
                continue
            end
            
            # locateArcAbove(event, y)
            middle = at(beachline, rel(event.site.point.x)...; mode = :linear) do arc
                # deleteEvent(arc)
                if !isnothing(arc.event)
                    delete!(events, arc.event)
                    arc.event = nothing
                end

                # breakArc(arc)
                a = Arc(arc.site)
                b = Arc(event.site)
                c = Arc(arc.site)

                a.left = arc.left
                c.right = arc.right
                
                return Replace(arc, b), InsertBefore(b, a), InsertAfter(b, c)
            end

            left = middle.prev
            right = middle.next

            # addEdge(left, middle)
            left.data.right = push!(voronoi, left.data.site.face)
            middle.data.left = push!(voronoi, middle.data.site.face)
            left.data.right.twin = middle.data.left
            middle.data.left.twin = left.data.right

            middle.data.right = middle.data.left
            right.data.left = left.data.right

            isnothing(left.prev) || add(left.prev.data, left.data, middle.data)
            isnothing(right.next) || add(middle.data, right.data, right.next.data)
        elseif event.type == :circle
            at(beachline, event.arc, -2:2, mode = :linear) do farleft, left, middle, right, farright
                # deleteEvent(left)
                if !isnothing(left.event)
                    delete!(events, left.event)
                    left.event = nothing
                end
                # deleteEvent(right)
                if !isnothing(right.event)
                    delete!(events, right.event)
                    right.event = nothing
                end

                # setDestination(left, middle)
                # setDestination(middle, right)
                vertex = 
                    left.right.origin = 
                    middle.left.destination = 
                    middle.right.origin = 
                    right.left.destination = 
                        push!(voronoi, event.point)

                middle.left.next = middle.right
                middle.right.prev = middle.left

                prev = left.right
                next = right.left

                # addEdge(left, right)
                left.right = push!(voronoi, left.site.face)
                right.left = push!(voronoi, right.site.face)
                left.right.twin = right.left
                right.left.twin = left.right

                # setOrigin(left, right)
                left.right.destination = right.left.origin = vertex

                # setPrevHalfEdge(left.right, prev)
                left.right.next = prev
                prev.prev = left.right
                # setPrevHalfEdge(next, right.left)
                next.next = right.left
                right.left.prev = next

                isnothing(farleft) || add(farleft, left, right)
                isnothing(farright) || add(left, right, farright)

                return Delete(middle)
            end
        else
            error("unknown event")
        end
    end

    for ver in voronoi.vertices
        expand!(bounds, ver.point)
    end

    links = Vector{LinkedVertex{T}}()
    cells = Dict{Int, Vector{Union{Nothing, LinkedVertex{T}}}}()
    LinkedRedBlackTree.rule(beachline, 0:1, mode = :linear) do left, right
        dir = ortho(left.site.point - right.site.point)
        origin = (left.site.point + right.site.point) / 2
        hit = raycast(bounds, origin, dir)

        # setDestination(left, right)
        vertex = 
            left.right.origin = 
            right.left.destination =
                push!(voronoi, hit.point)

        !haskey(cells, left.site.index) && (cells[left.site.index] = Vector{Union{Nothing, LinkedVertex{T}}}(nothing, 8))
        !haskey(cells, right.site.index) && (cells[right.site.index] = Vector{Union{Nothing, LinkedVertex{T}}}(nothing, 8))

        side = findfirst(==(hit.side), SIDES) - 1
        l = LinkedVertex(vertex, nothing, left.right)
        push!(links, l)
        cells[left.site.index][2side + 2] = l
        l = LinkedVertex(vertex, right.left, nothing)
        push!(links, l)
        cells[right.site.index][2side + 1] = l
    end

    for (_, cell) in cells
        for i in 0:4
            side = mod(i, 4)
            next = mod(i + 1, 4)
            if isnothing(cell[2side + 1]) && !isnothing(cell[2side + 2])
                prev = mod(i - 1, 4)
                l = LinkedVertex(push!(voronoi, bounds, SIDES[side + 1]), nothing, nothing)
                push!(links, l)
                cell[2prev + 2] = cell[2side + 1] = l
            elseif !isnothing(cell[2side + 1]) && isnothing(cell[2side + 2])
                prev = mod(i - 1, 4)
                l = LinkedVertex(push!(voronoi, bounds, SIDES[next + 1]), nothing, nothing)
                push!(links, l)
                cell[2side + 2] = cell[2next + 1] = l
            end
        end
    end

    for (i, cell) in cells
        for side in 0:3
            isnothing(cell[2side + 1]) && continue
            he = push!(voronoi, voronoi.faces[i])
            he.origin = cell[2side + 1].vertex
            he.destination = cell[2side + 2].vertex
            cell[2side + 1].next = he
            he.prev = cell[2side + 1].prev
            !isnothing(cell[2side + 1].prev) && (cell[2side + 1].prev.next = he)
            cell[2side + 2].prev = he
            he.next = cell[2side + 2].next
            !isnothing(cell[2side + 2].next) && (cell[2side + 2].next.prev = he)
        end
    end

    return voronoi, bounds
end

function bound!(voronoi::Voronoi{T}, bounds::AABB{T}) where T
    link = (a, a_side, b, b_side) -> begin
        he = a
        a_i = findfirst(==(a_side), SIDES) - 1
        b_i = findfirst(==(b_side), SIDES) - 1
        while a_i != b_i
            a_i = mod(a_i + 1, 4)
            he.next = push!(voronoi, he.parent)
            he.next.prev = he
            he.next.origin = he.destination
            he.next.destination = push!(voronoi, bounds, SIDES[a_i + 1])
            he = he.next
        end
        he.next = push!(voronoi, he.parent)
        he.next.prev = he
        b.prev = he.next
        he.next.next = b
        he.next.origin = he.destination
        he.next.destination = b.origin
    end

    closedlist = Set{HalfEdge{T, Face{T}}}()
    deletequeue = Vector{Vertex{T}}()
    for site in voronoi.sites
        he = site.face.children
        inside = isinside(bounds, he.origin.point)
        dirty = !inside
        incoming::Union{Nothing, HalfEdge{T, Face{T}}} = nothing
        outgoing::Union{Nothing, HalfEdge{T, Face{T}}} = nothing
        side_in = side_out = :none
        while true
            hits = raycastall(bounds, he.origin.point, he.destination.point)
            nextinside = isinside(bounds, he.destination.point)
            next = he.next
            if !inside && !nextinside
                push!(deletequeue, he.origin)
                if isempty(hits)
                    delete!(voronoi, he)
                elseif length(hits) == 2
                    if he.twin in closedlist
                        he.origin = he.twin.destination
                        he.destination = he.twin.origin
                    else
                        he.origin = push!(voronoi, hits[1].point)
                        he.destination = push!(voronoi, hits[2].point)
                    end
                    !isnothing(outgoing) && link(outgoing, side_out, incoming, side_in)
                    if isnothing(incoming)
                        incoming = he
                        side_in = hits[1].side
                    end
                    outgoing = he
                    side_out = hits[2].side
                else
                    error("this shoudn't happen")
                end
                push!(closedlist, he)
            elseif inside && !nextinside
                length(hits) == 1 || error("this shoudn't happen")
                if he.twin in closedlist
                    he.destination = he.twin.origin
                else
                    he.destination = push!(voronoi, hits[1].point)
                end
                outgoing = he
                side_out = hits[1].side
                push!(closedlist, he)
            elseif !inside && nextinside
                length(hits) == 1 || error("this shoudn't happen")
                push!(deletequeue, he.origin)
                if he.twin in closedlist
                    he.origin = he.twin.destination
                else
                    he.origin = push!(voronoi, hits[1].point)
                end
                !isnothing(outgoing) && link(outgoing, side_out, he, hits[1].side)
                if isnothing(incoming)
                    incoming = he
                    side_in = hits[1].side
                end
                push!(closedlist, he)
            end
            he = next
            inside = nextinside
            he === site.face.children && break
        end
        dirty && !isnothing(incoming) && link(outgoing, side_out, incoming, side_in)
        dirty && (site.face.children = incoming)
    end

    for ver in deletequeue
        delete!(voronoi, ver)
    end

    return nothing
end

function draw(v::Voronoi, bounds::AABB, filename; W = 500, H = 500, border = 5)
    Drawing(W + 2border, H + 2border, filename)
    
    scale = Point(W / (bounds.max.x - bounds.min.x), H / (bounds.max.y - bounds.min.y))
    offset1 = Point(bounds.min...)
    offset2 = Point(border, border)
    trafo = p -> (Point(p...) .- offset1) .* scale .+ offset2

    setcolor("black")
    box(trafo(bounds.min), trafo(bounds.max), action = :stroke)

    setcolor("black")
    for s in v.sites
        circle(trafo(s.point), 5, action = :fill)
    end

    setcolor("blue")
    for he in v.halfedges
        line(trafo(he.origin.point), trafo(he.destination.point), action = :stroke)
    end
        
    setcolor("red")
    for ver in v.vertices
        circle(trafo(ver.point), 5, action = :fill)
    end
    finish()
end

function main(N)
    points = [Vec2(rand(2)...) for i in 1:N]
    bounds = AABB{Float64}(-0.1, -0.1, 1.1, 1.1)
    voronoi, _ = fortune(points, bounds)
    bound!(voronoi, bounds)
    draw(voronoi, bounds, "fortune/voronoi.svg")
    
    #weights = rand(N)
    #voronoi, _ = fortune(collect(zip(points, weights)), bounds)
    #bound!(voronoi, bounds)
    #draw(voronoi, bounds, "fortune/weighted.svg")
end

main(20)