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
Flatten.flattenable(::Type{<:SFCCSolver}, ::Any) = false
"""
    SFCC{F,S} <: FreezeCurve

Generic representation of the soil freezing characteristic curve along with a nonlinear `solver`
for resolving the temperature-energy conservation law. The shape and parameters of the curve are
determined by the implementation of SFCCFunction `f`.
"""
struct SFCC{F,S} <: FreezeCurve
    f::F # freeze curve function f: (T,...) -> θ
    solver::S # solver for H -> T or T -> H
    SFCC(f::F, s::S=SFCCPreSolver()) where {F<:SFCCFunction,S<:SFCCSolver} = new{F,S}(f,s)
end
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
    temperature_residual(f::F, f_args::Fargs, hc, L, H, T) where {F<:SFCCFunction,Fargs<:Tuple}
    
Helper function for updating θw, C, and the residual. `hc` should be a function `θw -> C`
which computes the heat capacity from liquid water content (`θw`).
"""
@inline function temperature_residual(f::F, f_kwargs::NamedTuple, hc, L, H, T, residual_only=false) where {F<:SFCCFunction}
    θw = f(T; f_kwargs...)
    C = hc(θw)
    Tres = T - (H - θw*L) / C
    return residual_only ? Tres : (;Tres, θw, C)
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
    T;
    θtot=f.swrc.water.θtot,
    θsat=f.swrc.water.θsat,
    θres=f.swrc.water.θres, 
    Tₘ=f.freezethaw.Tₘ,
    α=f.swrc.α,
    n=f.swrc.n
)
    let θsat = max(θtot, θsat),
        g = f.g,
        Lsl = f.freezethaw.Lsl,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        ψ₀ = f.swrc(inv, θtot, θsat, θres, α, n),
        Tstar = Tₘ + g*Tₘ/Lsl*ψ₀,
        # matric potential as a function of T (Dall'Amico)
        ψ = ψ₀ + Lsl/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T);
        return f.swrc(ψ, θsat, θres, α, n)
    end
end
"""
    DallAmicoSalt{Tftp,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction

Freeze curve from Dall'Amico (2011) with saline freezing point depression.

Angelopoulos M, Westermann S, Overduin P, Faguet A, Olenchenko V, Grosse G, Grigoriev MN. Heat and salt flow in subsea permafrost
    modeled with CryoGRID2. Journal of Geophysical Research: Earth Surface. 2019 Apr;124(4):920-37.
"""
Base.@kwdef struct DallAmicoSalt{Tftp,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    saltconc::Tsc = 890.0u"mol/m^3" # salt concentration
    R::TR = 8.314459u"J/K/mol" #[J/K mol] universal gas constant
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
# DallAmico freeze curve with salt
function (f::DallAmicoSalt)(
    T;
    θtot=f.swrc.water.θtot,
    θsat=f.swrc.water.θsat,
    θres=f.swrc.water.θres, 
    Tₘ=f.freezethaw.Tₘ,
    saltconc=f.saltconc,
    α=f.swrc.α,
    n=f.swrc.n
)
    let θsat = max(θtot, θsat),
        g = f.g,
        R = f.R,
        ρw = f.swrc.water.ρw,
        Lsl = f.freezethaw.Lsl,
        Lf = Lsl*ρw,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        # freezing point depression based on salt concentration
        Tstar = Tₘ + Tₘ/Lf*(-R * saltconc * Tₘ),
        ψ₀ = f.swrc(inv, θtot, θsat, θres, α, n),
        # water matric potential
        ψ = ψ₀ + Lf / (ρw * g * Tstar) * (T - Tstar) * heaviside(Tstar-T),
        ψ = IfElse.ifelse(ψ < 0.0, ψ, 0.0);
        # van Genuchten evaulation to get θw
        return f.swrc(ψ,θres,θsat,α,n)
    end
end
# method dispatch to get SWRC for DallAmico freeze curves
swrc(f::Union{DallAmico,DallAmicoSalt}) = f.swrc
# use water properties from SWRC for DallAmico
SoilWaterProperties(f::Union{DallAmico,DallAmicoSalt}) = SoilWaterProperties(swrc(f))
"""
    McKenzie{Tftp,Twp,Tγ} <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
Base.@kwdef struct McKenzie{Tftp,Twp,Tγ} <: SFCCFunction
    freezethaw::Tftp = SoilFreezeThawProperties()
    water::Twp = SoilWaterProperties()
    γ::Tγ = 0.1u"K"
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
    δ::Tδ = 0.1u"K"
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
