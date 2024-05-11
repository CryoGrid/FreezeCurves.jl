abstract type SaturationState{T} end

struct MatricPotential{T} <: SaturationState{T}
    value::T
end

struct Saturation{T} <: SaturationState{T}
    value::T
end

"""
    SWRC

Base type for soil water retention curves (SWRC) which relate soil water matric potential to water content.
"""
abstract type SWRC end

"""
    Base.inv(f::SWRC)

Returns a function `(args...) -> f(inv, args...)` which, when implemented, evaluates the inverse of the
soil water retention curve.
"""
Base.inv(f::SWRC) = (args...; kwargs...) -> f(inv, args...; kwargs...)

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
    SoilWaterVolume(f::SWRC)

Get the `SoilWaterVolume` defined for the given `SWRC` `f`. Must be implemented
for all `SWRC` types. Default implementation retrieves the field `water`.
"""
SoilWaterVolume(f::SWRC) = f.vol

"""
    inflectionpoint(f::SWRC)

Returns the analytical solution for the inflection point (i.e. where ∂²θ/∂ψ² = 0) of the SWRC, if available.
"""
inflectionpoint(f::SWRC) = error("not implemented for $f")

"""
    waterpotential(f::SWRC, θ; θsat=f.vol.θsat, θres=f.vol.θres, kwargs...)

Returns the water potential at volumetric water content `θ` by evaluating the inverse of `f`.
"""
function waterpotential(f::SWRC, θ; θsat=f.vol.θsat, θres=f.vol.θres, kwargs...)
    let θsat = max(θ, θsat);
        f(inv, θ; θsat, θres, kwargs...)
    end
end

"""
    VanGenuchten{Tvol,Tα,Tn} <: SWRC

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct VanGenuchten{Tvol,Tα,Tn} <: SWRC
    vol::Tvol = SoilWaterVolume()
    α::Tα = 1.0u"1/m" # domain: (0,Inf)
    n::Tn = 2.0 # domain: (1,Inf)
end
function (f::VanGenuchten)(ψ; θsat=f.vol.θsat, θres=f.vol.θres, α=f.α, n=f.n)
    let m = 1-1/n,
        x = 1 + (-α*ψ)^n;
        # @assert isnan(x) || x > zero(x) "van Genuchten violated constraint $x > 0 with input $ψ and parameters θsat=$θsat, θres=$θres, α=$α, n=$n"
        IfElse.ifelse(ψ <= zero(ψ), θres + (θsat - θres)*x^(-m), θsat*one(ψ))
    end
end
function (f::VanGenuchten)(
    ::typeof(inv),
    θ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    α=f.α,
    n=f.n
)
    let m = 1-1/n;
        IfElse.ifelse(θ < θsat, -1/α*(((θ-θres)/(θsat-θres))^(-1/m)-1.0)^(1/n), zero(1/α)*θ)
    end
end

inflectionpoint(f::VanGenuchten; α=f.α, n=f.n) = -1/α*((n - 1) / n)^(1/n)

# Preset Van Genucthen curves from soil types;
# based on Carsel and Parrish (1988) but with maximum 5% residual water content
VanGenuchten(name::String) = VanGenuchten(Symbol(lowercase(replace(name, " " => ""))))
VanGenuchten(name::Symbol) = VanGenuchten(Val{name}())
VanGenuchten(::Val{:sand}) = VanGenuchten(α=0.074u"cm^-1", n=2.96, vol=SoilWaterVolume(θsat=0.33, θres=0.03))
VanGenuchten(::Val{:sandyloam}) = VanGenuchten(α=0.075u"cm^-1", n=2.28, vol=SoilWaterVolume(θsat=0.41, θres=0.05))
VanGenuchten(::Val{:siltloam}) = VanGenuchten(α=0.02u"cm^-1", n=1.41, vol=SoilWaterVolume(θsat=0.45, θres=0.05))
VanGenuchten(::Val{:sandyclayloam}) = VanGenuchten(α=0.06u"cm^-1", n=1.48, vol=SoilWaterVolume(θsat=0.39, θres=0.05))
VanGenuchten(::Val{:silt}) = VanGenuchten(α=0.016u"cm^-1", n=1.37, vol=SoilWaterVolume(θsat=0.46, θres=0.03))
VanGenuchten(::Val{:sandyclay}) = VanGenuchten(α=0.027u"cm^-1", n=1.23, vol=SoilWaterVolume(θsat=0.39, θres=0.05))
VanGenuchten(::Val{:siltyclay}) = VanGenuchten(α=0.005u"cm^-1", n=1.09, vol=SoilWaterVolume(θsat=0.38, θres=0.05))
VanGenuchten(::Val{:clay}) = VanGenuchten(α=0.008u"cm^-1", n=1.09, vol=SoilWaterVolume(θsat=0.38, θres=0.05))
VanGenuchten(::Val{:organic}) = VanGenuchten(α=0.013u"cm^-1", n=1.20, vol=SoilWaterVolume(θsat=0.77, θres=0.01))

"""
    BrooksCorey{Twp,Tψₛ,Tλ} <: SWRC

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct BrooksCorey{Tvol,Tψₛ,Tλ} <: SWRC
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
    IfElse.ifelse(ψ < -ψₛ, θres + (θsat - θres)*(-ψₛ / ψ)^λ, θsat*one(ψ))
end
function (f::BrooksCorey)(
    ::typeof(inv),
    θ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    ψₛ=f.ψₛ,
    λ=f.λ,
)
    IfElse.ifelse(θ < θsat, -ψₛ*((θ - θres)/(θsat - θres))^(-1/λ), -ψₛ*one(θ))
end
inflectionpoint(f::BrooksCorey; ψₛ=f.ψₛ) = ψₛ
