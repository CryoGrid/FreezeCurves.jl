"""
    SWRCFunction

Base type for soil water retention curves (SWRC) which relate soil water matric potential to water content.
"""
abstract type SWRCFunction end
"""
    Base.inv(f::SWRCFunction)

Returns a function `(args...) -> f(inv, args...)` which, when implemented, evaluates the inverse of the
soil water retention curve.
"""
Base.inv(f::SWRCFunction) = (args...; kwargs...) -> f(inv, args...; kwargs...)
"""
    SoilWaterProperties{Tρw,Tθres,Tθtot,Tθsat}

Struct containing basic physical constants and propeties related to soil pore water.
"""
Base.@kwdef struct SoilWaterProperties{Tρw,Tθres,Tθtot,Tθsat}
    ρw::Tρw = 1000.0u"kg/m^3" # density of water
    θres::Tθres = 0.0 # residual water content
    θtot::Tθtot = 0.5 # total water content
    θsat::Tθsat = 0.5 # saturated water content
    function SoilWaterProperties(ρw, θres, θtot, θsat)
        @assert θres < θtot <= θsat <= one(θsat)
        new{typeof(ρw),typeof(θres),typeof(θtot),typeof(θsat)}(ρw, θres, θtot, θsat)
    end
end
"""
    SoilWaterProperties(f::SWRCFunction)

Get the `SoilWaterProperties` defined for the given `SWRCFunction` `f`. Must be implemented
for all `SWRCFunction` types. Default implementation retrieves the field `water`.
"""
SoilWaterProperties(f::SWRCFunction) = f.water
"""
    waterpotential(f::SWRCFunction, θ; θsat=f.water.θsat, θres=f.water.θres, kwargs...)

Returns the water potential at volumetric water content `θ` by evaluating the inverse of `f`.
"""
function waterpotential(f::SWRCFunction, θ; θsat=f.water.θsat, θres=f.water.θres, kwargs...)
    let θsat = max(θ, θsat);
        f(inv, θ; θsat, θres, kwargs...)
    end
end
"""
    VanGenuchten{Twp,Tα,Tn} <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct VanGenuchten{Twp,Tα,Tn} <: SWRCFunction
    water::Twp = SoilWaterProperties()
    α::Tα = Param(1.0, units=u"1/m", domain=OpenInterval(0,Inf))
    n::Tn = Param(2.0, domain=Interval{:closed,:open}(1,Inf))
end
function (f::VanGenuchten)(ψ; θsat=stripparams(f.water.θsat), θres=stripparams(f.water.θres), α=stripparams(f.α), n=stripparams(f.n))
    let m = 1-1/n;
        IfElse.ifelse(ψ <= zero(ψ), θres + (θsat - θres)*(1 + (-α*ψ)^n)^(-m), θsat)
    end
end
function (f::VanGenuchten)(::typeof(inv), θ; θsat=stripparams(f.water.θsat), θres=stripparams(f.water.θres), α=stripparams(f.α), n=stripparams(f.n))
    let m = 1-1/n;
        IfElse.ifelse(θ < θsat, -1/α*(((θ-θres)/(θsat-θres))^(-1/m)-1.0)^(1/n), zero(1/α))
    end
end
"""
    BrooksCorey{Twp,Tψₛ,Tλ} <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct BrooksCorey{Twp,Tψₛ,Tλ} <: SWRCFunction
    water::Twp = SoilWaterProperties()
    ψₛ::Tψₛ = Param(0.01, units=u"m")
    λ::Tλ = Param(0.2, domain=OpenInterval(0,Inf))
end
function (f::BrooksCorey)(ψ; θsat=stripparams(f.water.θsat), θres=stripparams(f.water.θres), ψₛ=stripparams(f.ψₛ), λ=stripparams(f.λ))
    IfElse.ifelse(ψ < -ψₛ, θres + (θsat - θres)*(-ψₛ / ψ)^λ, θsat)
end
function (f::BrooksCorey)(::typeof(inv), θ; θsat=stripparams(f.water.θsat), θres=stripparams(f.water.θres), ψₛ=stripparams(f.ψₛ), λ=stripparams(f.λ))
    IfElse.ifelse(θ < θsat, -ψₛ*((θ - θres)/(θsat - θres))^(-1/λ), -ψₛ)
end
