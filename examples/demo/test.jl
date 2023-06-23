msg = "Hello World"
println(msg)

using Plots
f(x) = x^2
g(x) = x^3

plot(f)
plot(g)

α = 2
s = "aa"
split(s)
strip(s)
join([s,s],"s")
replace(s,"toSearch"=>"toReplace")

d = 'a'
fg = string(s,s)
a = "$s alice"

a = zeros(5)
b = [1,2,3]
T = Int64
a = T[]
x = [10, "foo", false]

c = zeros(2,3)
c[1:2, :]
a = [[1,2,3] [2,3,4]]
size(a)

a = (1,2,3)
aNamedTuple = (a=1, b=2)

myDict = Dict('a'=>1, 'b'=>2, 'c'=>3)
myDict['a']

for (k,v) in myDict
    println("$k is $v")
end

# mutable struct MyOwnType
#     property1
#     property2::String
# end

mutable struct MyOwnType{T<:Number}
    property1
    property2::String
    property3::T
end


abstract type MyOwnGenericAbstractType end
abstract type MyOwnAbstractType <: MyOwnGenericAbstractType end
mutable struct AConcreteType <: MyOwnAbstractType
  property1
  property2::String
end

myObject = MyOwnType("something","something",10)
a = myObject.property3 # 10


