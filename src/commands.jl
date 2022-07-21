export rule, at
export Insert, Delete, Replace, InsertBefore, InsertAfter

struct Command{ArgsType <: Tuple}
    type::Symbol
    args::ArgsType

    Command(type::Symbol, args...) = new{typeof(args)}(type, args)
end
Insert(el...) = Command(:insert, el...)
Delete(el...) = Command(:delete, el...)
Replace(old, new) = Command(:replace, old, new)
InsertBefore(pivot, el...) = Command(:insert_before, pivot, el...)
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

function at(callback::Function, tree::LRBTree, data, dis::AbstractVector{Int} = [0]; mode = :binary)
    node = search(tree, data; mode)
    isnothing(node) && return nothing

    execute!(callback(tree[node, dis]...), tree; mode)

    return node
end

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

function rule(callback::Function, tree::LRBTree, dis::AbstractVector{Int} = [0]; skip_border = true, mode = :binary)
    isempty(tree) && return nothing

    lo = skip_border ? max(1 - minimum(dis), 1) : 1
    hi = skip_border ? min(length(tree) - maximum(dis), length(tree)) : length(tree)
    hi < lo && return nothing

    cb = CommandBuffer()
    node = tree.min
    for i in 1:hi
        i >= lo && push!(cb, callback(tree[node, dis]...))
        node = node.next
    end
    execute!(cb, tree; mode)
end