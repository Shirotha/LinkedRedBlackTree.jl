export rule, at
export Insert, Delete, Replace, InsertBefore, InsertAfter

struct Command{ArgsType <: Tuple}
    type::Symbol
    args::ArgsType

    Command(type::Symbol, args...) = new{typeof(args)}(type, args)
end
"Insert new elements into a `LRBTree`"
Insert(el...) = Command(:insert, el...)
"Delete existing elements from a `LRBTree`"
Delete(el...) = Command(:delete, el...)
"""
Replace the data of an existing `LRBTree` element.

WARNING: This won't check if the order is correct!
Consider accessing the tree using the `mode = :linear` option for `at` and `rule`.
"""
Replace(old, new) = Command(:replace, old, new)
"""
Insert new elements before an existing `LRBTree` element.

WARNING: This won't check if the order is correct!
Consider accessing the tree using the `mode = :linear` option for `at` and `rule`.
"""
InsertBefore(pivot, el...) = Command(:insert_before, pivot, el...)
"""
Insert new elements after an existing `LRBTree` element.

WARNING: This won't check if the order is correct!
Consider accessing the tree using the `mode = :linear` option for `at` and `rule`.
"""
InsertAfter(pivot, el...) = Command(:insert_after, pivot, el...)

struct CommandBuffer
    commands::Vector{Command}

    CommandBuffer() = new(Vector{Command}())
end
function CommandBuffer(cmd::Union{Command, Tuple{Vararg{Command}}, <:AbstractVector{Command}})
    cb = CommandBuffer()
    push!(cb, cmd)
    return cb
end

Base.push!(cb::CommandBuffer, cmd::Union{Command, Tuple{Vararg{Command}}, <:AbstractVector{Command}}) = !isnothing(cmd) && (cmd isa Command ? push!(cb.commands, cmd) : append!(cb.commands, cmd))
Base.push!(cb::CommandBuffer, cmd) = nothing

function execute!(cb::CommandBuffer, tree::LRBTree; mode = :binary)
    for cmd in cb.commands
        if cmd.type == :insert
            insert!(tree, cmd.args...)
        elseif cmd.type == :delete
            delete!(tree, cmd.args...; mode)
        elseif cmd.type == :replace
            node = search(tree, first(cmd.args); mode)
            node.data = last(cmd.args)
        elseif cmd.type == :insert_before
            node = search(tree, first(cmd.args); mode)
            for d in Iterators.tail(cmd.args)
                insert_before!(tree, node, d)
            end
        elseif cmd.type == :insert_after
            node = search(tree, first(cmd.args); mode)
            for d in Iterators.reverse(Iterators.tail(cmd.args))
                insert_after!(tree, node, d)
            end
        else
            error("unknown command $(cmd.type)")
        end
    end
    return tree
end
function execute!(cmd::Union{Command, Tuple{Vararg{Command}}, AbstractVector{Command}}, tree::LRBTree; mode = :binary)
    cb = CommandBuffer(cmd)
    execute!(cb, tree; mode)
end
"""
Search the `LRBTree` for given `data` and run the `callback` on elements relative to the found element.

# Parameters
- `callback`: will only be callend when `data` has been found, arguments to this will be in the same order as `dis`
- `tree`: the `LRBTree` to search through
- `data`: the data to search for
- `dis`: indices relative to the found element, used for arguments of `callback`
- `mode`: when set to `:binary` (default) the tree will be searched using binary search, use `:linear` to avoid this (e.g. after using `InsertBefore`/`InsertAfter`)

When the return value of `callback` is `Command` (e.g. `Insert`, `Delete`) or list of them then they will be executed after the callback.
"""
function at(callback::Function, tree::LRBTree, data, dis::AbstractVector{Int} = [0]; mode = :binary)
    node = search(tree, data; mode)
    isnothing(node) && return nothing

    execute!(callback(tree[node, dis]...), tree; mode)

    return node
end
"""
Search the `LRBTree` with custom comparators and run the `callback` on elements relative to the found element.

# Parameters
- `callback`: will only be callend when `data` has been found, arguments to this will be in the same order as `dis`
- `tree`: the `LRBTree` to search through
- `dis`: indices relative to the found element, used for arguments of `callback`
- `mode`: when set to `:binary` (default) the tree will be searched using binary search, use `:linear` to avoid this (e.g. after using `InsertBefore`/`InsertAfter`)

# Comparators
`lt` and `gt` take the data of two adjacent elements (left element always first) and should return a `Bool`.

The tree will be traversed depending on those return values:

| lt \\ gt | `true`   | `false` |
|:---------|:---------|:--------|
| `true`   | go left  | go left |
| `false`  | go right | found   |

When the search reaches the far left/right of the tree that element will be returned.

When the return value of `callback` is `Command` (e.g. `Insert`, `Delete`) or list of them then they will be executed after the callback.
"""
function at(callback::Function, tree::LRBTree, lt::Function, gt::Function, dis::AbstractVector{Int} = [0]; mode = :binary)
    node = tree.root
    while true
        if !isnothing(node.prev) && lt(node.prev.data, node.data)
            node = node.left
        elseif !isnothing(node.next) && gt(node.data, node.next.data)
            node = node.right
        else
            break
        end
    end
    execute!(callback(tree[node, dis]...), tree; mode)

    return node
end
"""
Iterate over the `LRBTree` and run the `callback` on elements relative to the current element.

# Parameters
- `callback`: arguments to this will be in the same order as `dis` (can be `nothing` when `skip_border` is `false`)
- `tree`: the `LRBTree` to iterate over
- `dis`: indices relative to the found element, used for arguments of `callback`
- `skip_border`: when iteration elements that would require access to invalid indices should be skipped (invalid elements will be `nothing` when this is set to `false`)
- `mode`: when set to `:binary` (default) the tree will be searched using binary search, use `:linear` to avoid this (e.g. after using `InsertBefore`/`InsertAfter`)

When the return value of `callback` is `Command` (e.g. `Insert`, `Delete`) or list of them then they will be executed after all callbacks are done.
"""
function rule(callback::Function, tree::LRBTree, dis::AbstractVector{Int} = [0]; skip_border = true, mode = :binary)
    isempty(tree) && return nothing

    lo = skip_border ? max(1 - minimum(dis), 1) : 1
    hi = skip_border ? min(length(tree) - maximum(dis), length(tree)) : length(tree)
    hi < lo && return 0

    cb = CommandBuffer()
    node = tree.min
    for i in 1:hi
        i >= lo && push!(cb, callback(tree[node, dis]...))
        node = node.next
    end
    execute!(cb, tree; mode)

    return hi - lo + 1
end