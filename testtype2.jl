__precompile__()

module TestTypes

using JLD: save

import Base: push!

export TestType1, TestType2, TestType3, TestType4

struct TestType1
    value::Int
end

struct TestType2
    data::Vector{TestType1}

    TestType2() = new(TestType1[])
end

struct TestType3
    data::Vector{TestType2}

    TestType3() = new(TestType2[])
end

struct TestType4
    data::Vector{TestType3}

    TestType4() = new(TestType3[])
end

for n in 2:4
    @eval function Base.push!(test_type::$(Symbol(:TestType, n)), test_arg::$(Symbol(:TestType, n-1)))
        push!(test_type.data, test_arg)
    end
end

function main()
    test1 = TestType1(2)
    test2 = TestType2()
    test3 = TestType3()
    test4 = TestType4()

    push!(test2, test1)
    push!(test3, test2)
    push!(test4, test3)

    file = "test_resutls.jld"
    save(file, "test4", test4)
    info("Test results saved in: $file")
    return @show test4
end

end
