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

Returns a function that evaluates the inverse of the soil water retention curve, `f`.
"""
Base.inv(f::SWRC) = (args...; kwargs...) -> f(inv, args...; kwargs...)

"""
    derivative(f::SWRC)

Returns a function that evaluates the derivative `∂f∂ψ` of the soil water retention curve, `f`.
"""
derivative(f::SWRC) = (args...; kwargs...) -> f(derivative, args..., kwargs...)

"""
    SoilWaterVolume{Tρw,Tθres,Tθsat}

Struct containing basic physical constants and propeties related to soil pore water.
"""
Base.@kwdef struct SoilWaterVolume{Tρw,Tθres,Tθsat}
    ρw::Tρw = 1000.0u"kg/m^3" # density of water
    θres::Tθres = 0.0 # residual water content
    θsat::Tθsat = 0.5 # saturated water content
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
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj1980.03615995004400050002x.
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
        ifelse(ψ < zero(ψ), θres + (θsat - θres)*x^(-m), θsat*one(ψ))
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
    # one(r) instead of a Float64 literal, so that the result has the number type of the arguments
    let m = 1-1/n,
        r = (θ-θres)/(θsat-θres);
        ifelse(θ < θsat, -1/α*(r^(-1/m)-one(r))^(1/n), zero(1/α)*θ)
    end
end

function (f::VanGenuchten)(
    ::typeof(derivative),
    ψ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    α=f.α,
    n=f.n
)
    # evaluate the unsaturated branch at a negative potential where saturated (ψ ≥ 0), so that
    # neither branch involves powers of zero, which give NaN derivatives (and adjoints) for n < 2
    let m = 1-1/n,
        unsaturated = ψ < zero(ψ),
        ψ_unsat = ifelse(unsaturated, ψ, -oneunit(ψ)),
        x = 1 + (-α*ψ_unsat)^n,
        ∂θ∂ψ = α*n*m*(θsat - θres)*x^(-m-1)*(-α*ψ_unsat)^(n-1);
        ifelse(unsaturated, ∂θ∂ψ, zero(∂θ∂ψ))
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

R. H. Brooks and A. T. Corey. 1964. Hydraulic Properties of Porous Media and Their Relation to Drainage Design. Transactions of the ASAE, 7(1): 0026–0028. DOI: 10.13031/2013.40684
"""
Base.@kwdef struct BrooksCorey{Tvol,Tψₛ,Tλ} <: SWRC
    vol::Tvol = SoilWaterVolume()
    ψₛ::Tψₛ = 0.01u"m" # air entry suction
    λ::Tλ = 0.25 # domain: (0,Inf)
end

function (f::BrooksCorey)(
    ψ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    ψₛ=f.ψₛ,
    λ=f.λ
)
    ifelse(ψ < -ψₛ, θres + (θsat - θres)*(-ψₛ / ψ)^λ, θsat*one(ψ))
end

function (f::BrooksCorey)(
    ::typeof(inv),
    θ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    ψₛ=f.ψₛ,
    λ=f.λ,
)
    ifelse(θ < θsat, -ψₛ*((θ - θres)/(θsat - θres))^(-1/λ), -ψₛ*one(θ))
end

function (f::BrooksCorey)(
    ::typeof(derivative),
    ψ;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    ψₛ=f.ψₛ,
    λ=f.λ
)
    # evaluate the unsaturated branch below the air entry potential -ψₛ where saturated, and return a zero
    # derivative of the same type there
    let unsaturated = ψ < -ψₛ,
        ψ_unsat = ifelse(unsaturated, ψ, oftype(ψ, -ψₛ - oneunit(ψ))),
        ∂θ∂ψ = (θsat - θres)*λ*(-ψₛ / ψ_unsat)^(λ-1)*(ψₛ*ψ_unsat^-2);
        ifelse(unsaturated, ∂θ∂ψ, zero(∂θ∂ψ))
    end
end

inflectionpoint(f::BrooksCorey; ψₛ=f.ψₛ) = ψₛ

# Preset Brooks-Corey curves from soil types;
# based on Table 2 of Clapp and Hornberger (1978) with λ = 1/b and no residual water content:
# Clapp RB and Hornberger GM, 1978. Empirical equations for some soil hydraulic properties.
#     Water Resources Research, 14(4): 601–604. DOI: 10.1029/WR014i004p00601.
BrooksCorey(name::String) = BrooksCorey(Symbol(lowercase(replace(name, " " => ""))))
BrooksCorey(name::Symbol) = BrooksCorey(Val{name}())
BrooksCorey(::Val{:sand}) = BrooksCorey(ψₛ=12.1u"cm", λ=1/4.05, vol=SoilWaterVolume(θsat=0.395))
BrooksCorey(::Val{:loamysand}) = BrooksCorey(ψₛ=9.0u"cm", λ=1/4.38, vol=SoilWaterVolume(θsat=0.410))
BrooksCorey(::Val{:sandyloam}) = BrooksCorey(ψₛ=21.8u"cm", λ=1/4.90, vol=SoilWaterVolume(θsat=0.435))
BrooksCorey(::Val{:siltloam}) = BrooksCorey(ψₛ=78.6u"cm", λ=1/5.30, vol=SoilWaterVolume(θsat=0.485))
BrooksCorey(::Val{:loam}) = BrooksCorey(ψₛ=47.8u"cm", λ=1/5.39, vol=SoilWaterVolume(θsat=0.451))
BrooksCorey(::Val{:sandyclayloam}) = BrooksCorey(ψₛ=29.9u"cm", λ=1/7.12, vol=SoilWaterVolume(θsat=0.420))
BrooksCorey(::Val{:siltyclayloam}) = BrooksCorey(ψₛ=35.6u"cm", λ=1/7.75, vol=SoilWaterVolume(θsat=0.477))
BrooksCorey(::Val{:sandyclay}) = BrooksCorey(ψₛ=15.3u"cm", λ=1/10.4, vol=SoilWaterVolume(θsat=0.426))
BrooksCorey(::Val{:siltyclay}) = BrooksCorey(ψₛ=49.0u"cm", λ=1/10.4, vol=SoilWaterVolume(θsat=0.492))
BrooksCorey(::Val{:clay}) = BrooksCorey(ψₛ=40.5u"cm", λ=1/11.4, vol=SoilWaterVolume(θsat=0.482))
