module InterfaceMechanics

export hello, domath

hello(who::String) = "Hello, $who"
domath(x::Number) = x + 5

println("ladataan paketti")
greet() = print("Hello World!")

end # module
