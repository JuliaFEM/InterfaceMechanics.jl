# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/InterfaceMechanics.jl/blob/master/LICENSE
module InterfaceMechanics

using LinearAlgebra, Parameters, StaticArrays, ForwardDiff

abstract type AbstractInterface end
abstract type AbstractInterfaceState end

@generated function Base.:+(state::T, dstate::T) where {T <: AbstractInterfaceState}
   expr = [:(state.$p+ dstate.$p) for p in fieldnames(T)]
   return :(T($(expr...)))
end

function update!(interface::I) where {I <: AbstractInterface}
    interface.drivers += interface.ddrivers
    interface.parameters += interface.dparameters
    interface.variables = interface.variables_new
    reset!(interface)
end

function reset!(interface::I) where {I <: AbstractInterface}
    interface.ddrivers = typeof(interface.ddrivers)()
    interface.dparameters = typeof(interface.dparameters)()
    interface.variables_new = typeof(interface.variables_new)()
end

function integrate_interface!(interface::I) where {I<:AbstractInterface}
    error("One needs to define how to integrate interface $I!")
end

export integrate_interface!

export update!, reset!

include("couloumb.jl")
export Couloumb, CouloumbDriverState, CouloumbParameterState, CouloumbVariableState
# export hello, domath
#
# hello(who::String) = "Hello, $who"
# domath(x::Number) = x + 5
#
# println("ladataan paketti")
# greet() = print("Hello World!")

end # module
