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
    SoilWaterProperties(f::SFCCFunction)

Retrieves the nested `SoilWaterProperties` from the `SoilFreezeThawProperties` of the freeze curve `f`.
"""
SoilWaterProperties(f::SFCCFunction) = f.water
"""
    temperature_residual(f::F, f_args::Fargs, hc, L, H, T, ::Val{return_all}=Val{true}()) where {F<:SFCCFunction,Fargs<:Tuple,return_all}
    
Helper function for updating θw, C, and the residual. `hc` should be a function `θw -> C`
which computes the heat capacity from liquid water content (`θw`).
"""
@inline function temperature_residual(f::F, f_kwargs::NamedTuple, hc, L, H, T, ::Val{return_all}=Val{true}()) where {F<:SFCCFunction,return_all}
    θw = f(T; f_kwargs...)
    C = hc(θw)
    Tres = T - (H - θw*L) / C
    return return_all ? (;Tres, θw, C) : Tres
end
"""
    PainterKarra{Tftp,Tβ,Tω,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Painter SL, Karra S. Constitutive Model for Unfrozen Water Content in Subfreezing Unsaturated Soils. Vadose zone j. 2014 Apr;13(4):1-8.
"""
Base.@kwdef struct PainterKarra{Tftp,Tβ,Tω,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    β::Tβ = Param(1.0, domain=OpenInterval(0,Inf))
    ω::Tω = Param(1/β, domain=0..(1/β))
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
@inline function (f::PainterKarra)(
    T,
    ψ₀=nothing,
    ::Val{return_all}=Val{false}();
    θtot=f.swrc.water.θtot,
    θsat=f.swrc.water.θsat,
    θres=f.swrc.water.θres, 
    Tₘ=f.freezethaw.Tₘ,
    β=f.β,
    ω=f.ω,
    swrc_kwargs...
) where return_all
    _waterpotential(ψ₀, θtot; kwargs...) = ψ₀
    _waterpotential(::Nothing, θtot; kwargs...) = waterpotential(f.swrc, θtot; kwargs...)
    let θsat = max(θtot, θsat),
        ψ₀ = _waterpotential(ψ₀, θtot; θsat, θres, swrc_kwargs...),
        g = f.g,
        Lsl = f.freezethaw.Lsl,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        Tstar = Tₘ + ω*g*Tₘ/Lsl*ψ₀,
        # matric potential as a function of T (same as Dall'Amico but with β parameter)
        ψ = ψ₀ + β*Lsl/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T),
        θw = f.swrc(ψ; θsat, θres, swrc_kwargs...);
        if return_all
            return (; θw, ψ, Tstar)
        else
            return θw
        end
    end
end
"""
    DallAmico{Tftp,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Dall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.
"""
Base.@kwdef struct DallAmico{Tftp,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
@inline function (f::DallAmico)(
    T,
    ψ₀=nothing,
    ::Val{return_all}=Val{false}();
    θtot=f.swrc.water.θtot,
    θsat=f.swrc.water.θsat,
    θres=f.swrc.water.θres, 
    Tₘ=f.freezethaw.Tₘ,
    swrc_kwargs...
) where return_all
    # Dall'Amico is a special case of Painter-Karra with ω = β = 1
    pkfc = PainterKarra()
    ω = 1.0
    β = 1.0
    return pkfc(T, ψ₀, Val{return_all}(); θtot, θsat, θres, Tₘ, ω, β, swrc_kwargs...)
end
"""
    DallAmicoSalt{Tftp,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Freeze curve from Dall'Amico (2011) with saline freezing point depression.

Angelopoulos M, Westermann S, Overduin P, Faguet A, Olenchenko V, Grosse G, Grigoriev MN. Heat and salt flow in subsea permafrost
    modeled with CryoGRID2. Journal of Geophysical Research: Earth Surface. 2019 Apr;124(4):920-37.
"""
Base.@kwdef struct DallAmicoSalt{Tftp,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    saltconc::Tsc = Param(890.0, units=u"mol/m^3", domain=Interval{:closed,:open}(0,Inf)) # salt concentration
    R::TR = 8.314459u"J/K/mol" #[J/K mol] universal gas constant
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
# DallAmico freeze curve with salt
function (f::DallAmicoSalt)(
    T,
    ψ₀=nothing,
    ::Val{return_all}=Val{false}();
    θtot=f.swrc.water.θtot,
    θsat=f.swrc.water.θsat,
    θres=f.swrc.water.θres, 
    Tₘ=f.freezethaw.Tₘ,
    saltconc=f.saltconc,
    α=f.swrc.α,
    n=f.swrc.n
) where return_all
    _waterpotential(ψ₀, θtot; kwargs...) = ψ₀
    _waterpotential(::Nothing, θtot; kwargs...) = waterpotential(f.swrc, θtot; kwargs...)
    let θsat = max(θtot, θsat),
        ψ₀ = _waterpotential(ψ₀, θtot; θsat, θres, α, n),
        g = f.g,
        R = f.R,
        ρw = f.swrc.water.ρw,
        Lsl = f.freezethaw.Lsl,
        Lf = Lsl*ρw,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        # freezing point depression based on salt concentration
        Tstar = Tₘ + Tₘ/Lf*(-R * saltconc * Tₘ),
        # water matric potential
        ψ = ψ₀ + Lf / (ρw * g * Tstar) * (T - Tstar) * heaviside(Tstar-T),
        ψ = IfElse.ifelse(ψ < zero(ψ), ψ, zero(ψ));
        # van Genuchten evaulation to get θw
        θw = f.swrc(ψ; θres, θsat, α, n)
        if return_all
            return (; θw, ψ, Tstar)
        else
            return θw
        end
    end
end
# method dispatches for SWRC-based freeze curves
swrc(f::Union{DallAmico,DallAmicoSalt,PainterKarra}) = f.swrc
SoilWaterProperties(f::Union{DallAmico,DallAmicoSalt,PainterKarra}) = SoilWaterProperties(swrc(f))
"""
    McKenzie{Tftp,Twp,Tγ} <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
Base.@kwdef struct McKenzie{Tftp,Twp,Tγ} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    water::Twp = SoilWaterProperties()
    γ::Tγ = Param(0.1, units=u"K", domain=OpenInterval(0,Inf))
end
function (f::McKenzie)(
    T;
    θtot=f.water.θtot,
    θsat=f.water.θsat,
    θres=f.water.θres,
    Tₘ=f.freezethaw.Tₘ,
    γ=f.γ
)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θsat = max(θtot, θsat);
        return IfElse.ifelse(T <= Tₘ, θres + (θsat-θres)*exp(-((T-Tₘ)/γ)^2), θtot)
    end
end
"""
    Westermann{Tftp,Twp,Tδ} <: SFCCFunction

Westermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of
    wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,
    https://doi.org/10.5194/tc-5-945-2011, 2011. 
"""
Base.@kwdef struct Westermann{Tftp,Twp,Tδ} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    water::Twp = SoilWaterProperties()
    δ::Tδ = Param(0.1, units=u"K", domain=OpenInterval(0,Inf))
end
function (f::Westermann)(
    T;
    θtot=f.water.θtot,
    θsat=f.water.θsat,
    θres=f.water.θres,
    Tₘ=f.freezethaw.Tₘ,
    δ=f.δ
)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θsat = max(θtot, θsat);
        return IfElse.ifelse(T <= Tₘ, θres - (θsat-θres)*(δ/(T-Tₘ-δ)), θtot)
    end
end
