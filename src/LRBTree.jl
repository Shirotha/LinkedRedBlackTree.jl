export LRBTree, value, indexof

"""
Node in `LRBTree`.
In addition to the tree structure (`parent`, `left`, `right`) this is also doubly linked (`prev`, `next`).
"""
mutable struct LRBNode{T}
    data::Union{Nothing, T}
    "true = red, false = black"
    color::Bool

    parent::Union{Nothing, LRBNode{T}}
    left::Union{Nothing, LRBNode{T}}
    right::Union{Nothing, LRBNode{T}}

    prev::Union{Nothing, LRBNode{T}}
    next::Union{Nothing, LRBNode{T}}

    LRBNode{T}() where T = new{T}(nothing, true, nothing, nothing, nothing, nothing, nothing)
    LRBNode(d::T) where T = new{T}(d, true, nothing, nothing, nothing, nothing, nothing)
end
"""
Returns the data held by a `LRBNode`.
"""
value(node::LRBNode) = node.data
value(::Nothing) = nothing
(node::LRBNode)() = node.data

"""
Red Black Tree data structure that is also doubly linked.
This allows for duplicate keys (delete! will remove the first matching key).
"""
# NOTE: inheriting from AbstractArray adds a lot of unnessesary copying for some reason
mutable struct LRBTree{T, LT <: Function, EQ <: Function}# <: AbstractVector{T}
    nil::LRBNode{T}

    root::LRBNode{T}

    min::LRBNode{T}
    max::LRBNode{T}    

    count::Int

    lt::LT
    eq::EQ

    function LRBTree{T}(lt::LT, eq::EQ) where {T, LT <: Function, EQ <: Function}
        nil = LRBNode{T}()
        nil.color = false
        return new{T, LT, EQ}(nil, nil, nil, nil, 0, lt, eq)
    end
end
"""
Create a new `LRBTree`.

# Parameters
- `lt`: comparator function used to traverse the tree while performing binary search, arguments are the result of `by`
- `eq`: comparator function used to determinte equality for searching and deleting elements, arguments are the elements themselves
- `by`: this will be called on the elements before being passed to `lt`
- `rev`: reverses the order of the `LRBTree` (`false` = smallest element first)
"""
LRBTree{T}(; lt = isless, eq = (==), by = identity, rev::Bool = false) where T =
    LRBTree{T}(@inline((a, b) -> xor(rev, lt(by(a), by(b)))), eq)

function rotate_left!(tree::LRBTree{T}, node_x::LRBNode{T}) where T
    node_y = node_x.right
    node_x.right = node_y.left
    node_y.left !== tree.nil && (node_y.left.parent = node_x)
    node_y.parent = node_x.parent
    if isnothing(node_x.parent)
        tree.root = node_y
    elseif node_x === node_x.parent.left
        node_x.parent.left = node_y
    else
        node_x.parent.right = node_y
    end
    node_y.left = node_x
    node_x.parent = node_y
    return nothing
end

function rotate_right!(tree::LRBTree{T}, node_x::LRBNode{T}) where T
    node_y = node_x.left
    node_x.left = node_y.right
    node_y.right !== tree.nil && (node_y.right.parent = node_x)
    node_y.parent = node_x.parent
    if isnothing(node_x.parent)
        tree.root = node_y
    elseif node_x === node_x.parent.left
        node_x.parent.left = node_y
    else
        node_x.parent.right = node_y
    end
    node_y.right = node_x
    node_x.parent = node_y
    return nothing
end

function fix_insert!(tree::LRBTree{T}, node::LRBNode{T}) where T
    parent::Union{Nothing, LRBNode{T}} = nothing
    grand_parent::Union{Nothing, LRBNode{T}} = nothing
    while node !== tree.root && node.parent.color
        parent = node.parent
        grand_parent = parent.parent
        isnothing(grand_parent) && return nothing
        if parent === grand_parent.left
            uncle = grand_parent.right
            if uncle.color
                grand_parent.color = true
                parent.color = false
                uncle.color = false
                node = grand_parent
            else
                node === parent.right && (node = parent; rotate_left!(tree, node))
                parent = node.parent
                grand_parent = parent.parent
                parent.color = false
                grand_parent.color = true
                rotate_right!(tree, grand_parent)
            end
        else
            uncle = grand_parent.left
            if uncle.color
                grand_parent.color = true
                parent.color = false
                uncle.color = false
                node = grand_parent
            else
                node === parent.left && (node = parent; rotate_right!(tree, node))
                parent = node.parent
                grand_parent = parent.parent
                parent.color = false
                grand_parent.color = true
                rotate_left!(tree, grand_parent)
            end
        end
    end
    tree.root.color = false
    return nothing
end

function fix_delete!(tree::LRBTree{T}, node::Union{Nothing, LRBNode{T}}) where T
    while node !== tree.root && !node.color
        if node === node.parent.left
            sibling = node.parent.right
            if sibling.color
                sibling.color = false
                node.parent.color = true
                rotate_left!(tree, node.parent)
                sibling = node.parent.right
            end
            if !sibling.right.color && !sibling.left.color
                sibling.color = true
                node = node.parent
            else
                if !sibling.right.color
                    sibling.left.color = false
                    sibling.color = true
                    rotate_right!(tree, sibling)
                    sibling = node.parent.right
                end
                sibling.color = node.parent.color
                node.parent.color = false
                sibling.right.color = false
                rotate_left!(tree, node.parent)
                node = tree.root
            end
        else
            sibling = node.parent.left
            if sibling.color
                sibling.color = false
                node.parent.color = true
                rotate_right!(tree, node.parent)
                sibling = node.parent.left
            end
            if !sibling.right.color && !sibling.left.color
                sibling.color = true
                node = node.parent
            else
                if !sibling.left.color
                    sibling.right.color = false
                    sibling.color = true
                    rotate_left!(tree, sibling)
                    sibling = node.parent.left
                end
                sibling.color = node.parent.color
                node.parent.color = false
                sibling.left.color = false
                rotate_right!(tree, node.parent)
                node = tree.root
            end
        end
    end
    node.color = false
    return nothing
end

function insert_root!(tree::LRBTree{T}, node::LRBNode{T}) where T
    node.parent = nothing
    tree.root = node

    tree.min = node
    tree.max = node

    node.left = tree.nil
    node.right = tree.nil

    tree.count = 1

    fix_insert!(tree, node)

    return tree
end

function insert_left!(tree::LRBTree{T}, parent::LRBNode{T}, node::LRBNode{T}) where T
    node.parent = parent
    parent.left = node

    node.prev = parent.prev
    node.next = parent
    if isnothing(parent.prev) 
        tree.min = node
    else
        parent.prev.next = node
    end
    parent.prev = node

    node.left = tree.nil
    node.right = tree.nil

    tree.count += 1

    fix_insert!(tree, node)

    return tree
end

function insert_right!(tree::LRBTree{T}, parent::LRBNode{T}, node::LRBNode{T}) where T
    node.parent = parent
    parent.right = node

    node.prev = parent
    node.next = parent.next
    if isnothing(parent.next)
        tree.max = node
    else
        parent.next.prev = node
    end
    parent.next = node

    node.left = tree.nil
    node.right = tree.nil

    tree.count += 1

    fix_insert!(tree, node)

    return tree
end

function insert_before!(tree::LRBTree{T}, parent::LRBNode{T}, node::LRBNode{T}) where T
    if parent.left === tree.nil
        insert_left!(tree, parent, node)
    else
        insert_right!(tree, parent.prev, node)
    end

    return tree
end
insert_before!(tree::LRBTree{T}, parent::LRBNode{T}, data::T) where T = insert_before!(tree, parent, LRBNode(data))

function insert_after!(tree::LRBTree{T}, parent::LRBNode{T}, node::LRBNode{T}) where T
    if parent.right === tree.nil
        insert_right!(tree, parent, node)
    else
        insert_left!(tree, parent.next, node)
    end

    return tree
end
insert_after!(tree::LRBTree{T}, parent::LRBNode{T}, data::T) where T = insert_after!(tree, parent, LRBNode(data))

function insert_at!(tree::LRBTree{T}, parent::Union{Nothing, LRBNode{T}}, node::LRBNode{T}) where T
    if isnothing(parent)
          insert_root!(tree, node)
    elseif tree.lt(node.data, parent.data)
        insert_left!(tree, parent, node)
    else
        insert_right!(tree, parent, node)
    end

    return tree
end

"""
Deletes given `LRBNode` from a `LRBTree`.
"""
function Base.delete!(tree::LRBTree{T}, z::LRBNode{T}) where T

    function replace!(tree::LRBTree{T}, old::LRBNode{T}, new::LRBNode{T})
        if isnothing(old.parent)
            tree.root = new
        elseif old === old.parent.left
            old.parent.left = new
        else
            old.parent.right = new
        end
        new.parent = old.parent
        return nothing
    end

    function minimum(tree::LRBTree{T}, node::LRBNode{T})
        node === tree.nil && return node
        while node.left !== tree.nil
            node = node.left
        end
        return node
    end

    z === tree.nil && return tree

    if isnothing(z.prev)
        if isnothing(z.next)
            tree.min = tree.nil
            tree.max = tree.nil
        else
            tree.min = z.next
            z.next.prev = nothing
        end
    elseif isnothing(z.next)
        tree.max = z.prev
        z.prev.next = nothing
    else
        z.prev.next = z.next
        z.next.prev = z.prev
    end

    y = z
    y_color = y.color
    if z.left === tree.nil
        x = z.right
        replace!(tree, z, z.right)
    elseif z.right === tree.nil
        x = z.left
        replace!(tree, z, z.left)
    else
        y = minimum(tree, z.right)
        y_color = y.color
        x = y.right
        if y.parent === z
            x.parent = y
        else
            replace!(tree, y, y.right)
            y.right = z.right
            y.right.parent = y
        end

        replace!(tree, z, y)
        y.left = z.left
        y.left.parent = y
        y.color = z.color
    end

    y_color || fix_delete!(tree, x)

    tree.count -= 1

    return tree
end

"""
Searches for a `LRBNode` in a given `LRBTree`.
Equality is determined by `LRBTree.eq`.

# mode
- :binary - uses binary search
- :linear - avoids using `LRBTree.lt` and searches the `LRBTree` like an unordered list

"""
function search(tree::LRBTree{T}, data::T; mode = :binary) where T
    if mode == :binary
        node = tree.root
        
        while node !== tree.nil && !tree.eq(data, node.data)
            node = tree.lt(data, node.data) ? node.left : node.right
        end

        return node
    elseif mode == :linear
        node = tree.min
        while !isnothing(node) && !tree.eq(data, node.data)
            node = node.next
        end

        return node
    else
        error("unknown mode")
    end
end

"""
Inserts a `LRBNode` to the correct position of a `LRBTree`.
"""
function Base.insert!(tree::LRBTree{T}, node::LRBNode{T}) where T
    node_x = tree.root
    node_y::Union{Nothing, LRBNode{T}} = nothing
    while node_x !== tree.nil
        node_y = node_x
        node_x = tree.lt(node.data, node_x.data) ? node_x.left : node_x.right
    end

    insert_at!(tree, node_y, node)
    return tree
end

Base.insert!(tree::LRBTree{T}, data::T) where T = insert!(tree, LRBNode(data))
Base.insert!(tree::LRBTree{T}, data) where T = insert!(tree, convert(T, data))
Base.insert!(tree::LRBTree, data...) = foreach(d -> insert!(tree, d), data)

Base.push!(tree::LRBTree, data...) = insert!(tree, data...)
Base.append!(tree::LRBTree, data::AbstractVector...) = insert!(tree, vcat(data...)...)

function Base.delete!(tree::LRBTree{T}, data::T; mode = :binary) where T
    node = search(tree, data; mode)
    tree.eq(data, node.data) && delete!(tree, node)
    return tree
end    
Base.delete!(tree::LRBTree{T}, data; mode = :binary) where T = delete!(tree, convert(T, data); mode)
Base.delete!(tree::LRBTree, data...; mode = :binary) = foreach(d -> delete!(tree, d; mode), data)

Base.in(data::T, tree::LRBTree{T}) where T = tree.eq(data, search(tree, data).data)
Base.in(data, tree::LRBTree{T}) where T = applicable(convert, T, data) ? in(convert(T, data), tree) : false

Base.IndexStyle(::Type{<:LRBTree}) = IndexLinear()
Base.size(tree::LRBTree) = (tree.count,)
Base.length(tree::LRBTree) = tree.count
Base.axes(tree::LRBTree) = (Base.OneTo(tree.count),)
Base.firstindex(::LRBTree) = 1
Base.lastindex(tree::LRBTree) = tree.count
Base.checkindex(tree::LRBTree, i::Int) = 1 <= i <= tree.count
Base.checkbounds(tree::LRBTree, i) = checkindex(tree, i) || throw(BoundsError(tree, i))
function Base.getindex(tree::LRBTree, i::Int)
    @boundscheck checkbounds(tree, i)
    if i > tree.count ÷ 2
        for (j, el) in Iterators.reverse(enumerate(tree))
            i == j && return el
        end
    else
        for (j, el) in enumerate(tree)
            i == j && return el
        end
    end
end
Base.checkindex(tree::LRBTree, is::AbstractVector{Int}) = all(i -> checkindex(tree, i), is)
"""
Accesses multiple indices as a batch, only traversing the linked list once.
"""
function Base.getindex(tree::LRBTree{T}, is::AbstractVector{Int}) where T
    isempty(is) && return Vector{T}()

    @boundscheck checkbounds(tree, is)

    order = sortperm(is)
    sorted = @inbounds is[order]

    result = Vector{T}(undef, length(is))

    if last(sorted) > length(tree) - first(sorted) + 1
        i = length(is)
        for (el, j) in Iterators.reverse(enumerate(tree))
            if @inbounds sorted[i] == j
                @inbounds result[i] = el
                i -= 1
                i < 1 && break
            end
        end
    else
        i = 1
        for (el, j) in enumerate(tree)
            if @inbounds sorted[i] == j
                @inbounds result[i] = el
                i += 1
                i > length(is) && break
            end
        end
    end

    return @inbounds result[invperm(order)]
end
"""
Accesses multiple indices relative to a given `LRBNode`, invalid indices will cause `nothing` to be returned.
"""
function Base.getindex(tree::LRBTree{T}, pivot::LRBNode{T}, dis::AbstractVector{Int}) where T
    isempty(dis) && return Vector{Union{Nothing, LRBNode{T}}}()

    min = minimum(dis)
    max = maximum(dis)
    order = sortperm(dis, by = di -> di < 0 ? min - di : di)
    sorted = @inbounds dis[order]

    result = Vector{Union{Nothing, T}}(undef, length(dis))

    i = 1
    node = pivot
    for j in -1:-1:min
        !isnothing(node) && (node = node.prev)
        if @inbounds sorted[i] == j
            @inbounds result[i] = isnothing(node) ? nothing : node.data
            i += 1
        end
    end
    node = pivot
    for j in 0:max
        if @inbounds sorted[i] == j
            @inbounds result[i] = isnothing(node) ? nothing : node.data
            i += 1
        end
        !isnothing(node) && (node = node.next)
    end
    
    return @inbounds result[invperm(order)]
end

Base.first(tree::LRBTree) = tree.min.data
Base.last(tree::LRBTree) = tree.max.data
Base.minimum(tree::LRBTree) = tree.min.data
Base.maximum(tree::LRBTree) = tree.max.data

Base.iterate(tree::LRBTree) = tree.count == 0 ? nothing : (tree.min.data, tree.min)
Base.iterate(tree::LRBTree, current::LRBNode) = current === tree.max ? nothing : (current.next.data, current.next)
Base.iterate(reverse::Iterators.Reverse{<:LRBTree}) = let tree = reverse.itr
    tree.count == 0 ? nothing : (tree.max.data, tree.max)
end
Base.iterate(reverse::Iterators.Reverse{<:LRBTree}, current::LRBNode) = let tree = reverse.itr
    current === tree.min ? nothing : (current.prev.data, current.prev)
end

function Base.vcat(trees::LRBTree{T}...) where T
    n = 0
    for t in trees
        n += t.count
    end
    result = Vector{T}(undef, n)
    i = 1
    for t in trees, el in t
        @inbounds result[i] = el
        i += 1
    end
    return result
end
Base.collect(tree::LRBTree) = vcat(tree)

function LRBTree(tree::LRBTree{T}) where T
    copy = LRBTree{T}(tree.lt, tree.eq)
    append!(copy, tree)
    return copy
end
function LRBTree(vec::AbstractVector{T}; options...) where T
    tree = LRBTree{T}(; options...)
    append!(tree, vec)
    return tree
end
function LRBTree{T}(vec::AbstractVector; options...) where T
    tree = LRBTree{T}(; options...)
    append!(tree, vec)
    return tree
end
LRBTree(data...; options...) where T = LRBTree(collect(data); options...)

Base.convert(type::Type{LRBTree{T}}, tree::LRBTree) where T = LRBTree{T}(tree)
Base.convert(type::Type{<:LRBTree}, vec::AbstractVector) = LRBTree(vec)

function Base.pop!(tree::LRBTree)
    node = tree.min
    delete!(tree, node)
    return node.data
end
"""
Finds the index of a `LRBNode`.
"""
function indexof(node::LRBNode)
    i = 1
    current = node.prev
    while !isnothing(current)
        i += 1
        current = current.prev
    end
    return i
end
indexof(tree::LRBTree, data) = indexof(search(tree, data))

# TODO: custom broatcasting so that getindex is not called every time

function Base.showarg(io::IO, tree::LRBTree{T}, toplevel::Bool) where T
    toplevel && print(io, tree.count, "-element ")
    print(io, "LRBTree{", T, "}")
end
function Base.show(io::IO, tree::LRBTree{T}) where T
    print(io, "LRBTree{", T, "}(", collect(tree), ")")
end
function Base.show(io::IO, ::MIME"text/plain", tree::LRBTree{T}) where T
    if get(io, :compact, false)
        show(io, tree)
    else
        Base.showarg(io, tree, true)
        println(io, ":")
        Base.print_array(io, collect(tree))
        println(io)
    end
end