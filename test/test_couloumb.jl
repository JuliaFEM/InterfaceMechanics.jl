# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/InterfaceMechanics.jl/blob/master/LICENSE
using InterfaceMechanics, Test, StaticArrays

parameters = CouloumbParameterState(mu=0.2,
                                    elastic_slip=0.001,
                                    c0 = 0.001,
                                    p0 = 100.0)

h0 = -0.005 # Initial overclosure
dh = 0.006 # Overclosure
du = 1.0
drivers = CouloumbDriverState(displacements = SVector{3,Float64}([h0, 0.0, 0.0]))
dtime = 1.0
ddrivers = CouloumbDriverState(time = dtime,
                               displacements = SVector{3,Float64}([dh, 0.0, 0.0]))

interface = Couloumb(parameters=parameters,
                     drivers=drivers,
                     ddrivers=ddrivers)

@info("********* Starting simulation *********")
# Step 1: Initialize contact
@info "time = $(interface.drivers.time), traction = $(interface.variables.traction)"
integrate_interface!(interface)
update!(interface)
@info "time = $(interface.drivers.time), traction = $(interface.variables.traction)"
@test interface.variables.traction[1]>0
@test interface.variables.traction[2]==0
@test interface.variables.traction[3]==0

# Step 2: Tangential movement 1
ddrivers = CouloumbDriverState(time = dtime,
                               displacements = SVector{3,Float64}([0.0, du, 0.0]))
interface.ddrivers = ddrivers
integrate_interface!(interface)
update!(interface)
@info "time = $(interface.drivers.time), traction = $(interface.variables.traction)"
@test interface.variables.traction[1]>0
@test interface.variables.traction[2]>0
@test interface.variables.traction[3]==0

# Step 3: Tangential movement 2
ddrivers = CouloumbDriverState(time = dtime,
                               displacements = SVector{3,Float64}([0.0, -du, 0.0]))
interface.ddrivers = ddrivers
integrate_interface!(interface)
update!(interface)
@info "time = $(interface.drivers.time), traction = $(interface.variables.traction)"
@test interface.variables.traction[1]>0
@test interface.variables.traction[2]<0
@test interface.variables.traction[3]==0

# Step 3: Tangential movement 3
ddrivers = CouloumbDriverState(time = dtime,
                               displacements = SVector{3,Float64}([0.0, du/sqrt(2), du/sqrt(2)]))
interface.ddrivers = ddrivers
integrate_interface!(interface)
update!(interface)
@info "time = $(interface.drivers.time), traction = $(interface.variables.traction)"
@test interface.variables.traction[1]>0
@test interface.variables.traction[2]>0
@test interface.variables.traction[3]>0
