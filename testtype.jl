

type TestType1
    value::Int
end

type TestType2
    testvec::Any
    TestType2()=new(Test2[])
end

type TestType3
    testtype2::Any
    TestType3()=new(TestType2[])
end

type TestType4
    testtype3::Any
    TestType4()=new(TestType3[])
end


test1=Test2(2);
test2=TestType2();
test3=TestType3();
test4=TestType4();
push!(test2.testvec,test1)
push!(test3.testtype2,test2)
push!(test4.testtype3,test3)

using JLD
save("/Users/pmxct2/Documents/FosterDayanMorris/testype.jld", "test4",test4)
