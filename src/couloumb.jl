# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/InterfaceMechanics.jl/blob/master/LICENSE

@with_kw mutable struct CouloumbDriverState <: AbstractInterfaceState
    time :: Float64 = zero(Float64)
    displacements :: SVector{3,Float64} = zero(SVector{3,Float64})
end

@with_kw mutable struct CouloumbParameterState <: AbstractInterfaceState
    mu :: Float64 = 0.0
    elastic_slip :: Float64 = 0.000
    p0 :: Float64 = 0.0
    c0 :: Float64 = 0.0
end

@with_kw mutable struct CouloumbVariableState <: AbstractInterfaceState
    traction :: SVector{3,Float64} = zero(SVector{3,Float64})
    slip :: SVector{2,Float64} = zero(SVector{2,Float64})
    frictional_energy :: Float64 = 0.0
    jacobian :: SMatrix{3,3,Float64} = zero(SMatrix{3,3,Float64})
end

@with_kw mutable struct Couloumb <: AbstractInterface
    drivers :: CouloumbDriverState = CouloumbDriverState()
    ddrivers :: CouloumbDriverState = CouloumbDriverState()
    variables :: CouloumbVariableState = CouloumbVariableState()
    variables_new :: CouloumbVariableState = CouloumbVariableState()
    parameters :: CouloumbParameterState = CouloumbParameterState()
    dparameters :: CouloumbParameterState = CouloumbParameterState()
end

function integrate_interface!(interface::Couloumb)
    par = interface.parameters
    var = interface.variables
    ddri = interface.ddrivers
    dri = interface.drivers
    @unpack mu, elastic_slip, p0, c0 = par
    @unpack displacements, time = dri
    h0, u10, u20 = displacements
    dtime = ddri.time

    # Drivers at start state
    @unpack traction, slip, frictional_energy, jacobian = var

    # Begin integration
    function dvars(ddisplacement)
        dh, du1, du2 = ddisplacement
        h = h0 + dh
        # Define exponential pressure-overclosure
        p_func(h) = p0/(exp(1.0)-1.0)*((h/c0 + 1.0)*(exp(h/c0 + 1.0) - 1.0))
        dp_func(h) = ForwardDiff.derivative(p_func, h)
        if h <= -c0 # No contact
            p = 0.0
            tau1 = 0.0
            tau2 = 0.0
        elseif h<=6*c0 # Exponential pressure-overclosure
            p = p_func(h)
        else # Linearized continuation after h>6*c0
            p = p_func(6.0*c0) + dp_func(6.0*c0)*(h - 6.0*c0)
        end
        dp = p - traction[1]
        # Trial values for interface tangential movement
        u1 = u10 + du1
        u2 = u20 + du2
        tau_cr = mu*p
        gamma1_tr = u1 - slip[1]
        gamma2_tr = u2 - slip[2]
        gamma_eq_tr = sqrt(gamma1_tr^2 + gamma2_tr^2)
        k = tau_cr/elastic_slip
        if gamma_eq_tr == 0 # Avoid division by zero
            n1 = 1.0/sqrt(2.0)
            n2 = 1.0/sqrt(2.0)
        else
            n1 = gamma1_tr/gamma_eq_tr
            n2 = gamma2_tr/gamma_eq_tr
        end
        jacobian_ = zeros(3,3)
        if gamma_eq_tr>elastic_slip # Slipping
            dslip1 = n1*(gamma_eq_tr - elastic_slip)
            dslip2 = n2*(gamma_eq_tr - elastic_slip)
            tau1 = n1*tau_cr
            tau2 = n2*tau_cr
            dtau1 = tau1 - traction[2]
            dtau2 = tau2 - traction[3]
            dfrictional_energy_ = tau_cr*(gamma_eq_tr - elastic_slip)
            jacobian_[1,1] = dp_func(h)
            jacobian_[2,1] = n1*mu*jacobian_[1,1]
            jacobian_[2,2] = (1.0-n1*n1)*tau_cr/gamma_eq_tr
            jacobian_[2,3] = -n1*n2*tau_cr/gamma_eq_tr
            jacobian_[3,1] = n2*mu*jacobian_[1,1]
            jacobian_[3,2] = -n1*n2*tau_cr/gamma_eq_tr
            jacobian_[3,3] = (1.0-n2*n2)*tau_cr/gamma_eq_tr
        else # Stick/reversal
            dtau1 = k*du1
            dtau2 = k*du2
            dslip1 = 0.0
            dslip2 = 0.0
            dfrictional_energy_ = 0.0
            jacobian_[1,1] = dp_func(h)
            jacobian_[2,1] = (traction[2]+dtau1)/tau_cr*mu
            jacobian_[2,2] = k*jacobian_[1,1]
            jacobian_[3,1] = (traction[3]+dtau2)/tau_cr*mu
            jacobian_[3,3] = k*jacobian_[1,1]
        end
        return [dp, dtau1, dtau2], [dslip1, dslip2], dfrictional_energy_, jacobian_
    end
    dtraction, dslip, dfrictional_energy, jacobian = dvars(ddri.displacements)
    @info("ddisplacements = $(ddri.displacements)")
    @info("traction = $traction, slip = $slip, frictional_energy = $frictional_energy")
    @info("dtraction = $dtraction, dslip = $dslip, dfrictional_energy = $dfrictional_energy")
    #jacobian = ForwardDiff.jacobian(x -> dvars(x)[1], ddri.displacements)
    @info("jacobian = $jacobian")
    frictional_energy += dfrictional_energy
    traction += SVector{3, Float64}(dtraction)
    slip += SVector{2, Float64}(dslip)
    variables_new = CouloumbVariableState(traction=SVector{3, Float64}(traction),
                slip = SVector{2, Float64}(slip),
                frictional_energy = frictional_energy,
                jacobian = SMatrix{3,3,Float64}(jacobian))
    interface.variables_new = variables_new
end
