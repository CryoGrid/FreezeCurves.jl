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
derivative(f::SWRC) = (args...; kwargs...) -> f(derivative, args...; kwargs...)

"""
    inflectionpoint(f::SWRC)

Returns the analytical solution for the inflection point (i.e. where ∂²θ/∂ψ² = 0) of the SWRC, if available.
"""
inflectionpoint(f::SWRC) = error("not implemented for $f")

"""
    waterpotential(f::SWRC, θ; θsat=0.5, θres=0.0, kwargs...)

Returns the water potential at volumetric water content `θ` by evaluating the inverse of `f`.
"""
function waterpotential(f::SWRC, θ; θsat=0.5, θres=0.0, kwargs...)
    let θsat = max(θ, θsat);
        f(inv, θ; θsat, θres, kwargs...)
    end
end

"""
    VanGenuchten{Tα,Tn} <: SWRC

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj1980.03615995004400050002x.
"""
Base.@kwdef struct VanGenuchten{Tα,Tn} <: SWRC
    α::Tα = 1.0u"1/m" # domain: (0,Inf)
    n::Tn = 2.0 # domain: (1,Inf)
end
function (f::VanGenuchten)(ψ; θsat=0.5, θres=0.0, α=f.α, n=f.n)
    let m = 1-1/n,
        x = 1 + (-α*ψ)^n;
        # @assert isnan(x) || x > zero(x) "van Genuchten violated constraint $x > 0 with input $ψ and parameters θsat=$θsat, θres=$θres, α=$α, n=$n"
        ifelse(ψ < zero(ψ), θres + (θsat - θres)*x^(-m), θsat*one(ψ))
    end
end

function (f::VanGenuchten)(
    ::typeof(inv),
    θ;
    θsat=0.5,
    θres=0.0,
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
    θsat=0.5,
    θres=0.0,
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
VanGenuchten(::Val{:sand}) = VanGenuchten(α=0.074u"cm^-1", n=2.96)
VanGenuchten(::Val{:sandyloam}) = VanGenuchten(α=0.075u"cm^-1", n=2.28)
VanGenuchten(::Val{:siltloam}) = VanGenuchten(α=0.02u"cm^-1", n=1.41)
VanGenuchten(::Val{:sandyclayloam}) = VanGenuchten(α=0.06u"cm^-1", n=1.48)
VanGenuchten(::Val{:silt}) = VanGenuchten(α=0.016u"cm^-1", n=1.37)
VanGenuchten(::Val{:sandyclay}) = VanGenuchten(α=0.027u"cm^-1", n=1.23)
VanGenuchten(::Val{:siltyclay}) = VanGenuchten(α=0.005u"cm^-1", n=1.09)
VanGenuchten(::Val{:clay}) = VanGenuchten(α=0.008u"cm^-1", n=1.09)
VanGenuchten(::Val{:organic}) = VanGenuchten(α=0.013u"cm^-1", n=1.20)

"""
    BrooksCorey{Tψₛ,Tλ} <: SWRC

R. H. Brooks and A. T. Corey. 1964. Hydraulic Properties of Porous Media and Their Relation to Drainage Design. Transactions of the ASAE, 7(1): 0026–0028. DOI: 10.13031/2013.40684
"""
Base.@kwdef struct BrooksCorey{Tψₛ,Tλ} <: SWRC
    ψₛ::Tψₛ = 0.01u"m" # air entry suction
    λ::Tλ = 0.25 # domain: (0,Inf)
end

function (f::BrooksCorey)(
    ψ;
    θsat=0.5,
    θres=0.0,
    ψₛ=f.ψₛ,
    λ=f.λ
)
    ifelse(ψ < -ψₛ, θres + (θsat - θres)*(-ψₛ / ψ)^λ, θsat*one(ψ))
end

function (f::BrooksCorey)(
    ::typeof(inv),
    θ;
    θsat=0.5,
    θres=0.0,
    ψₛ=f.ψₛ,
    λ=f.λ,
)
    ifelse(θ < θsat, -ψₛ*((θ - θres)/(θsat - θres))^(-1/λ), -ψₛ*one(θ))
end

function (f::BrooksCorey)(
    ::typeof(derivative),
    ψ;
    θsat=0.5,
    θres=0.0,
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
BrooksCorey(::Val{:sand}) = BrooksCorey(ψₛ=12.1u"cm", λ=1/4.05)
BrooksCorey(::Val{:loamysand}) = BrooksCorey(ψₛ=9.0u"cm", λ=1/4.38)
BrooksCorey(::Val{:sandyloam}) = BrooksCorey(ψₛ=21.8u"cm", λ=1/4.90)
BrooksCorey(::Val{:siltloam}) = BrooksCorey(ψₛ=78.6u"cm", λ=1/5.30)
BrooksCorey(::Val{:loam}) = BrooksCorey(ψₛ=47.8u"cm", λ=1/5.39)
BrooksCorey(::Val{:sandyclayloam}) = BrooksCorey(ψₛ=29.9u"cm", λ=1/7.12)
BrooksCorey(::Val{:siltyclayloam}) = BrooksCorey(ψₛ=35.6u"cm", λ=1/7.75)
BrooksCorey(::Val{:sandyclay}) = BrooksCorey(ψₛ=15.3u"cm", λ=1/10.4)
BrooksCorey(::Val{:siltyclay}) = BrooksCorey(ψₛ=49.0u"cm", λ=1/10.4)
BrooksCorey(::Val{:clay}) = BrooksCorey(ψₛ=40.5u"cm", λ=1/11.4)

"""
    watercontent(::Type{NF}=Float64; θsat=NF(0.5), θres=NF(0.0)) where {NF}
    watercontent(texture::Symbol, source::Symbol)
    watercontent(texture::Val, source::Val=Val(:clapphornberger))

Returns a `NamedTuple` `(; θsat, θres)` of the saturated and residual volumetric water contents.

The first method constructs the `NamedTuple` directly from the keyword arguments, with the default values
`θsat = 0.5` and `θres = 0.0` converted to the number type `NF`.

The remaining methods look up the values for the given soil `texture` (e.g. `:sand`, `:siltloam`) from
the data `source`, given either as `Symbol`s or as `Val`s. Available sources:
- `:clapphornberger` (default for the `Val` method): Table 2 of Clapp and Hornberger (1978) with no
  residual water content; consistent with the `BrooksCorey` presets.
- `:carselparrish`: Carsel and Parrish (1988) with a maximum 5% residual water content;
  consistent with the `VanGenuchten` presets.

The result can be passed as keyword arguments when evaluating a curve, e.g.
`BrooksCorey(:sand)(ψ; watercontent(Val(:sand))...)` or
`VanGenuchten(:sand)(ψ; watercontent(:sand, :carselparrish)...)`.
"""
watercontent(::Type{NF} = Flaot64; θsat = NF(0.5), θres = NF(0.0)) where {NF} = (; θsat, θres)
watercontent(texture::Symbol, source::Symbol) = watercontent(Val(texture), Val(source))
watercontent(texture::Val) = watercontent(texture, Val(:clapphornberger))
# saturated and residual water contents from Carsel and Parrish (1988), consistent with the VanGenuchten presets
watercontent(::Val{:sand}, ::Val{:carselparrish}) = (θsat=0.33, θres=0.03)
watercontent(::Val{:sandyloam}, ::Val{:carselparrish}) = (θsat=0.41, θres=0.05)
watercontent(::Val{:siltloam}, ::Val{:carselparrish}) = (θsat=0.45, θres=0.05)
watercontent(::Val{:sandyclayloam}, ::Val{:carselparrish}) = (θsat=0.39, θres=0.05)
watercontent(::Val{:silt}, ::Val{:carselparrish}) = (θsat=0.46, θres=0.03)
watercontent(::Val{:sandyclay}, ::Val{:carselparrish}) = (θsat=0.39, θres=0.05)
watercontent(::Val{:siltyclay}, ::Val{:carselparrish}) = (θsat=0.38, θres=0.05)
watercontent(::Val{:clay}, ::Val{:carselparrish}) = (θsat=0.38, θres=0.05)
watercontent(::Val{:organic}, ::Val{:carselparrish}) = (θsat=0.77, θres=0.01)
# saturated and residual water contents from Clapp and Hornberger (1978), consistent with the BrooksCorey presets
watercontent(::Val{:sand}, ::Val{:clapphornberger}) = (θsat=0.395, θres=0.0)
watercontent(::Val{:loamysand}, ::Val{:clapphornberger}) = (θsat=0.410, θres=0.0)
watercontent(::Val{:sandyloam}, ::Val{:clapphornberger}) = (θsat=0.435, θres=0.0)
watercontent(::Val{:siltloam}, ::Val{:clapphornberger}) = (θsat=0.485, θres=0.0)
watercontent(::Val{:loam}, ::Val{:clapphornberger}) = (θsat=0.451, θres=0.0)
watercontent(::Val{:sandyclayloam}, ::Val{:clapphornberger}) = (θsat=0.420, θres=0.0)
watercontent(::Val{:siltyclayloam}, ::Val{:clapphornberger}) = (θsat=0.477, θres=0.0)
watercontent(::Val{:sandyclay}, ::Val{:clapphornberger}) = (θsat=0.426, θres=0.0)
watercontent(::Val{:siltyclay}, ::Val{:clapphornberger}) = (θsat=0.492, θres=0.0)
watercontent(::Val{:clay}, ::Val{:clapphornberger}) = (θsat=0.482, θres=0.0)
