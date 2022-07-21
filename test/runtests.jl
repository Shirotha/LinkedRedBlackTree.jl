using Test

include("../src/LinkedRedBlackTree.jl")
using .LinkedRedBlackTree

@testset verbose=true "LRBTree" begin
    @testset "Insertion" begin
        tree = LRBTree{Int}()
        insert!(tree, 2)
        insert!(tree, 1)
        insert!(tree, 3)
        @test length(tree) == 3
        @test tree.root.data == 2
        @test collect(tree) == [1, 2, 3] 
    end
    @testset "Deletion" begin
        tree = LRBTree(1:3)
        delete!(tree, 2)
        delete!(tree, 1)
        delete!(tree, 3)
        @test isempty(tree)
        @test tree.root == tree.nil
    end
    @testset "Search" begin
        tree = LRBTree(10:-2:2)
        @test 4 in tree
        @test !(3 in tree)
    end
    @testset "Same keys" begin
        tree = LRBTree(1, 1)
        @test collect(tree) == [1, 1]
        @test 1 in tree
        delete!(tree, 1)
        @test tree.count == 1 
    end
    @testset "Options" begin
        struct TestPoint
            x::Float64
            y::Float64
        end
        Base.isapprox(a::TestPoint, b::TestPoint) = isapprox(a.x, b.x) && isapprox(a.y, b.y)

        tree = LRBTree{TestPoint}(by = p -> p.y, eq = isapprox, rev = true)
        @test tree.lt(TestPoint(2, 2), TestPoint(4, 1))
        insert!(tree, TestPoint(1, 2), TestPoint(3, 0), TestPoint(-2, -2), TestPoint(0, 4))

        @test TestPoint(3, 0) in tree
        @test first(tree) ≈ TestPoint(0, 4)
        @test tree[end] ≈ TestPoint(-2, -2)
    end
    @testset "Indexing" begin
        tree = LRBTree(1:10)
        @test all(i -> tree[i] == i, 1:10)
        @test tree[[2, 7]] == [2, 7]
        @test tree[3:7] == 3:7
        @test tree[tree.min, -1:0] == [nothing, 1]
    end
    @testset "Commands" begin
        tree = LRBTree([1, 16])
        for _ in 1:5
            rule(tree, 0:1) do a, b
                if abs(b - a) > 1
                    Insert((a + b) ÷ 2)
                end
            end
        end
        @test all(tree .== 1:16)
    end
    @testset "Manual Ordering" begin
        tree = LRBTree([5])
        at(tree, 5) do e
            return Replace(e, 3), InsertBefore(3, 6), InsertAfter(3, 2)
        end
        @test collect(tree) == [6, 3, 2]
    end
end