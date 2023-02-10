"""
    SFCCFunction

Base type for a soil freeze characteristic curve (SFCC) function.
Subtypes should be callable structs that implement the freeze curve and contain
any necessary additional constants or configuration options. User-specified parameters
can either be supplied in the struct or declared as model parameters via the `variables`
method.
"""
abstract type SFCCFunction end
"""
    swrc(::SFCCFunction)

Retrieves the `SWRCFunction` embedded within this `SFCCFunction`, if defined.
The default value for freeze curve functions without an associated soil-water
retention curve is `nothing`.
"""
swrc(::SFCCFunction) = nothing
"""
    SFCCSolver

Base type representing non-linear solvers for implicit SFCC functions.
"""
abstract type SFCCSolver end
# never flatten SFCCSolver types
Flatten.flattenable(::Type{<:SFCCSolver}, ::Type{Val{fieldname}}) where {fieldname} = false
"""
    SFCCTable{F,I} <: SFCCFunction

Pre-tabulated soil freezing characteristic curve. See `Tabulated` and `tabulate`.
"""
struct SFCCTable{F,I} <: SFCCFunction
    f::F
    f_tab::I
end
(f::SFCCTable)(args...) = f.f_tab(args...)
"""
    Tabulated(f::SFCCFunction, args...)

Produces an `SFCCTable` function which is a tabulation of `f`.
"""
Tabulated(f::SFCCFunction, args...; kwargs...) = SFCCTable(f, tabulate(f, args...; kwargs...))
"""
    SoilFreezeThawProperties{TTₘ,TLsl}

Struct containing constants/parameters common to some or all SFCCs.
"""
Base.@kwdef struct SoilFreezeThawProperties{TTₘ,TLsl}
    Tₘ::TTₘ = 0.0u"°C" # melting temperature
    Lsl::TLsl = 3.34e5u"J/kg" # specific latent heat of fusion of water
end
"""
    SoilFreezeThawProperties(f::SFCCFunction)

Retrieves the default `SoilFreezeThawProperties` from `f`; should be defined for all freeze curves.
"""
SoilFreezeThawProperties(f::SFCCFunction) = f.freezethaw
"""
    SoilWaterVolume(f::SFCCFunction)

Retrieves the nested `SoilWaterVolume` from the `SoilFreezeThawProperties` of the freeze curve `f`.
"""
SoilWaterVolume(f::SFCCFunction) = f.vol
"""
    inflectionpoint(f::SFCCFunction)

Returns the analytical inflection point (i.e. where ∂²θ/∂T^2 = 0), if available.
"""
inflectionpoint(f::SFCCFunction) = error("not implemented for $f")
"""
    temperature_residual(f::F, f_args::Fargs, hc, L, H, T, ψ₀=nothing) where {F<:SFCCFunction}
    
Helper function for updating θw, C, and the residual. `hc` should be a function `θw -> C`
which computes the heat capacity from liquid water content (`θw`).
"""
@inline function temperature_residual(f::F, f_kwargs::NamedTuple, hc, L, H, T, sat=1.0) where {F<:SFCCFunction}
    # helper function to handle case where SFCC has no SWRC and thus does not accept ψ₀ as an argument
    _invoke_f(::Nothing, T, sat) = f(T; f_kwargs...)
    _invoke_f(::SWRCFunction, T, sat) = f(T, sat; f_kwargs...)
    θw = _invoke_f(swrc(f), T, sat)
    C = hc(θw)
    Tres = T - (H - θw*L) / C
    return Tres, θw, C
end
"""
    PainterKarra{TFT,Tβ,Tω,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Painter SL, Karra S. Constitutive Model for Unfrozen Water Content in Subfreezing Unsaturated Soils. Vadose zone j. 2014 Apr;13(4):1-8.
"""
Base.@kwdef struct PainterKarra{TFT,Tβ,Tω,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::TFT = SoilFreezeThawProperties()
    β::Tβ = 1.0 # domain: (0,Inf)
    ω::Tω = 1/β # domain: (0,1/β)
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
@inline function (f::PainterKarra)(
    T,
    sat=1.0,
    ::Val{output}=Val{:θw}();
    θsat=f.swrc.vol.θsat,
    θres=f.swrc.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    β=f.β,
    ω=f.ω,
    swrc_kwargs...
) where output
    let θtot = θres + (θsat-θres)*sat,
        ψ₀ = waterpotential(swrc(f), θtot; θsat, θres, swrc_kwargs...),
        g = f.g,
        β = β,
        ω = ω,
        Lsl = f.freezethaw.Lsl,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        Tstar = Tₘ + ω*g*Tₘ/Lsl*ψ₀;
        # matric potential as a function of T (same as Dall'Amico but with β parameter)
        ψ = ψ₀ + β*Lsl/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T)
        # output is a compile-time constant so will be compiled away
        if output == :all || output == true # allow 'true' for backwards compatibility w/ 0.4
            θw = f.swrc(ψ; θsat, θres, swrc_kwargs...)
            return (; θw, ψ, Tstar)
        elseif output == :θw || output == false # allow 'false' for backwards compatibility w/ 0.4
            θw = f.swrc(ψ; θsat, θres, swrc_kwargs...)
            return θw
        elseif output == :ψ
            return ψ
        elseif output == :Tstar
            return Tstar
        end
    end
end
function inflectionpoint(
    f::PainterKarra,
    sat=1.0;
    θsat=f.swrc.vol.θsat,
    θres=f.swrc.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    β=f.β,
    ω=f.ω,
    swrc_kwargs...
)
    let θtot = θres + (θsat-θres)*sat,
        ψstar = inflectionpoint(f.swrc; swrc_kwargs...),
        ψ₀ = waterpotential(swrc(f), θtot; θsat, θres, swrc_kwargs...),
        g = f.g,
        β = β,
        ω = ω,
        Lsl = f.freezethaw.Lsl,
        Tₘ = normalize_temperature(Tₘ),
        Tstar = Tₘ + ω*g*Tₘ/Lsl*ψ₀;
        return (ψstar - ψ₀)/(β*Lsl)*(g*Tstar) + Tstar
    end
end
"""
    DallAmico{TFT,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Dall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.
"""
Base.@kwdef struct DallAmico{TFT,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::TFT = SoilFreezeThawProperties()
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
@inline function (f::DallAmico)(
    T,
    sat=1.0,
    ::Val{output}=Val{:θw}();
    θsat=f.swrc.vol.θsat,
    θres=f.swrc.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    swrc_kwargs...
) where output
    # Dall'Amico is a special case of Painter-Karra with ω = β = 1
    pkfc = PainterKarra(freezethaw=f.freezethaw, g=f.g, swrc=f.swrc)
    ω = 1.0
    β = 1.0
    return pkfc(T, sat, Val{output}(); θsat, θres, Tₘ, ω, β, swrc_kwargs...)
end
function inflectionpoint(
    f::DallAmico,
    sat=1.0;
    θsat=f.swrc.vol.θsat,
    θres=f.swrc.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    swrc_kwargs...
)
    return inflectionpoint(PainterKarra(), sat; θsat, θres, Tₘ, swrc_kwargs...)
end
"""
    DallAmicoSalt{TFT,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Freeze curve from Dall'Amico (2011) with saline freezing point depression.

Angelopoulos M, Westermann S, Overduin P, Faguet A, Olenchenko V, Grosse G, Grigoriev MN. Heat and salt flow in subsea permafrost
    modeled with CryoGRID2. Journal of Geophysical Research: Earth Surface. 2019 Apr;124(4):920-37.
"""
Base.@kwdef struct DallAmicoSalt{TFT,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::TFT = SoilFreezeThawProperties()
    saltconc::Tsc = 890.0u"mol/m^3" # salt concentration
    R::TR = 8.314459u"J/K/mol" #[J/K mol] universal gas constant
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
# DallAmico freeze curve with salt
function (f::DallAmicoSalt)(
    T,
    sat=1.0,
    ::Val{output}=Val{:θw}();
    θsat=f.swrc.vol.θsat,
    θres=f.swrc.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    saltconc=f.saltconc,
    swrc_kwargs...
) where output
    let θtot = θres + (θsat-θres)*sat,
        ψ₀ = waterpotential(swrc(f), θtot; θsat, θres, swrc_kwargs...),
        g = f.g,
        R = f.R,
        ρw = f.swrc.vol.ρw,
        Lsl = f.freezethaw.Lsl,
        Lf = Lsl*ρw,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        # freezing point depression based on salt concentration
        Tstar = Tₘ + Tₘ/Lf*(-R * saltconc * Tₘ),
        # water matric potential
        ψ = ψ₀ + Lf / (ρw * g * Tstar) * (T - Tstar) * heaviside(Tstar-T),
        ψ = IfElse.ifelse(ψ < zero(ψ), ψ, zero(ψ));
        # output is a compile-time constant so will be compiled away
        if output == :all || output == true # allow 'true' for backwards compatibility w/ 0.4
            θw = f.swrc(ψ; θsat, θres, swrc_kwargs...)
            return (; θw, ψ, Tstar)
        elseif output == :θw || output == false # allow 'false' for backwards compatibility w/ 0.4
            θw = f.swrc(ψ; θsat, θres, swrc_kwargs...)
            return θw
        elseif output == :ψ
            return ψ
        elseif output == :Tstar
            return Tstar
        end
    end
end
# method dispatches for SWRC-based freeze curves
swrc(f::Union{DallAmico,DallAmicoSalt,PainterKarra}) = f.swrc
SoilWaterVolume(f::Union{DallAmico,DallAmicoSalt,PainterKarra}) = SoilWaterVolume(swrc(f))
"""
    McKenzie{TFT,Tvol,Tγ} <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
Base.@kwdef struct McKenzie{TFT,Tvol,Tγ} <: SFCCFunction
    freezethaw::TFT = SoilFreezeThawProperties()
    vol::Tvol = SoilWaterVolume()
    γ::Tγ = 0.1u"K" # domain (0,Inf)
end
function (f::McKenzie)(
    T,
    sat=1.0;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    γ=f.γ,
)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θtot = sat*θsat;
        return IfElse.ifelse(T <= Tₘ, θres + (θtot-θres)*exp(-((T-Tₘ)/γ)^2), θtot)
    end
end
"""
    Westermann{TFT,Tvol,Tδ} <: SFCCFunction

Westermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of
    wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,
    https://doi.org/10.5194/tc-5-945-2011, 2011. 
"""
Base.@kwdef struct Westermann{TFT,Tvol,Tδ} <: SFCCFunction
    freezethaw::TFT = SoilFreezeThawProperties()
    vol::Tvol = SoilWaterVolume()
    δ::Tδ = 0.1u"K" # domain: (0,Inf)
end
function (f::Westermann)(
    T,
    sat=1.0;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    δ=f.δ,
)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θtot = sat*θsat;
        return IfElse.ifelse(T <= Tₘ, θres - (θtot-θres)*(δ/(T-Tₘ-δ)), θtot)
    end
end

"""
    Hu2020{TFT,Tvol,Tb}  <: SFCCFunction

Soil freezing characteristic curve formulation of Hu et al. 2020.

Hu G, Zhao L, Zhu X, Wu X, Wu T, Li R, Xie C, Hao J. Review of algorithms and parameterizations to determine unfrozen water content in frozen soil. Geoderma. 2020 Jun 1;368:114277.
"""
Base.@kwdef struct Hu2020{TFT,Tvol,Tb} <: SFCCFunction
    freezethaw::TFT = SoilFreezeThawProperties()
    vol::Tvol = SoilWaterVolume()
    b::Tb = 0.01
end
function (f::Hu2020)(
    T,
    sat=1.0;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    Tₘ=f.freezethaw.Tₘ,
    b=f.b,
)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θtot = sat*θsat;
        return IfElse.ifelse(T <= Tₘ, θres + (θtot-θres)*(1 - ((Tₘ - T) / Tₘ)^b), θtot)
    end
end

"""
    PowerLaw{Tvol,Tα,Tβ} <: SFCCFunction

Commonly used power law relation of Lovell (1957).

Lovell, C.: Temperature effects on phase composition and strength of partially frozen soil, Highway Research Board Bulletin, 168, 74–95, 1957.
"""
Base.@kwdef struct PowerLaw{Tvol,Tα,Tβ} <: SFCCFunction
    vol::Tvol = SoilWaterVolume()
    α::Tα = 0.01
    β::Tβ = 0.5
end
function (f::PowerLaw)(
    T,
    sat=1.0;
    θsat=f.vol.θsat,
    θres=f.vol.θres,
    α=f.α,
    β=f.β,
)
    let T = ustrip(normalize_temperature(T)),
        Tₘ = ustrip(normalize_temperature(-α^(1/β))),
        θtot = sat*θsat;
        return IfElse.ifelse(T < Tₘ, θres + (θtot-θres)*(α*abs(T - Tₘ)^(-β)), θtot)
    end
end

const SUTRAIce_Exp = McKenzie
const SUTRAIce_Power = PowerLaw
