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
    SoilWaterVolume{Tρw,Tθres,Tθsat}

Struct containing basic physical constants and propeties related to soil pore water.
"""
Base.@kwdef struct SoilWaterVolume{Tρw,Tθres,Tθsat}
    ρw::Tρw = 1000.0u"kg/m^3" # density of water
    θres::Tθres = 0.0 # residual water content
    θsat::Tθsat = 0.5 # saturated water content
    function SoilWaterVolume(ρw, θres, θsat)
        @assert zero(θres) <= θres < θsat <= one(θsat)
        new{typeof(ρw),typeof(θres),typeof(θsat)}(ρw, θres, θsat)
    end
end
"""
    SoilWaterVolume(f::SWRCFunction)

Get the `SoilWaterVolume` defined for the given `SWRCFunction` `f`. Must be implemented
for all `SWRCFunction` types. Default implementation retrieves the field `water`.
"""
SoilWaterVolume(f::SWRCFunction) = f.vol
"""
    inflectionpoint(f::SWRCFunction)

Returns the analytical solution for the inflection point (i.e. where ∂²θ/∂ψ² = 0) of the SWRC, if available.
"""
inflectionpoint(f::SWRCFunction) = error("not implemented for $f")
"""
    waterpotential(f::SWRCFunction, θ; θsat=f.vol.θsat, θres=f.vol.θres, kwargs...)

Returns the water potential at volumetric water content `θ` by evaluating the inverse of `f`.
"""
function waterpotential(f::SWRCFunction, θ; θsat=f.vol.θsat, θres=f.vol.θres, kwargs...)
    let θsat = max(θ, θsat);
        f(inv, θ; θsat, θres, kwargs...)
    end
end
"""
    VanGenuchten{Tvol,Tα,Tn} <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct VanGenuchten{Tvol,Tα,Tn} <: SWRCFunction
    vol::Tvol = SoilWaterVolume()
    α::Tα = 1.0u"1/m" # domain: (0,Inf)
    n::Tn = 2.0 # domain: (1,Inf)
end
function (f::VanGenuchten)(ψ; θsat=f.vol.θsat, θres=f.vol.θres, α=f.α, n=f.n)
    let m = 1-1/n;
        IfElse.ifelse(ψ <= zero(ψ), θres + (θsat - θres)*(1 + (-α*ψ)^n)^(-m), θsat)
    end
end
function (f::VanGenuchten)(
    ::typeof(inv), θ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    α=f.α,
    n=f.n
)
    let m = 1-1/n;
        IfElse.ifelse(θ < θsat, -1/α*(((θ-θres)/(θsat-θres))^(-1/m)-1.0)^(1/n), zero(1/α))
    end
end
inflectionpoint(f::VanGenuchten; α=f.α, n=f.n) = -1/α*((n - 1) / n)^(1/n)
"""
    BrooksCorey{Twp,Tψₛ,Tλ} <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct BrooksCorey{Tvol,Tψₛ,Tλ} <: SWRCFunction
    vol::Tvol = SoilWaterVolume()
    ψₛ::Tψₛ = 0.01u"m"
    λ::Tλ = 0.2 # domain: (0,Inf)
end
function (f::BrooksCorey)(
    ψ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    ψₛ=f.ψₛ,
    λ=f.λ
)
    IfElse.ifelse(ψ < -ψₛ, θres + (θsat - θres)*(-ψₛ / ψ)^λ, θsat)
end
function (f::BrooksCorey)(
    ::typeof(inv), θ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    ψₛ=f.ψₛ,
    λ=f.λ,
)
    IfElse.ifelse(θ < θsat, -ψₛ*((θ - θres)/(θsat - θres))^(-1/λ), -ψₛ)
end
inflectionpoint(f::BrooksCorey; ψₛ=f.ψₛ) = ψₛ
