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
    SFCCSolver

Base type representing non-linear solvers for implicit SFCC functions.
"""
abstract type SFCCSolver end
"""
    SFCC{F,S} <: FreezeCurve

Generic representation of the soil freezing characteristic curve. The shape and parameters
of the curve are determined by the implementation of SFCCFunction `f`. For implicit functions,
a non-linear `solver` must also be specified.
"""
struct SFCC{F,S} <: FreezeCurve
    f::F # freeze curve function f: (T,...) -> θ
    solver::S # solver for H -> T or T -> H
    SFCC(f::F, s::S) where {F<:SFCCFunction,S<:Union{Nothing,SFCCSolver}} = new{F,S}(f,s)
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
    SoilFreezingProperties{TTₘ,TLsl,Tρw,Tθres,Tθtot,Tθsat}

Struct containing constants/parameters common to some or all SFCCs.
"""
Base.@kwdef struct SoilFreezingProperties{TTₘ,TLsl,Tρw,Tθres,Tθtot,Tθsat}
    Tₘ::TTₘ = 0.0u"°C" # melting temperature
    Lsl::TLsl = 3.34e5u"J/kg" # specific latent heat of fusion of water
    ρw::Tρw = 1000.0u"kg/m^3" # density of water
    θres::Tθres = 0.0 # residual water content
    θtot::Tθtot = 0.5 # total water content
    θsat::Tθsat = 0.5 # saturated water content
    function SoilFreezingProperties(Tₘ, Lsl, ρw, θres, θtot, θsat)
        @assert θres < θtot <= θsat <= one(θsat)
        new{typeof(Tₘ),typeof(Lsl),typeof(ρw),typeof(θres),typeof(θtot),typeof(θsat)}(Tₘ, Lsl, ρw, θres, θtot, θsat)
    end
end
"""
    temperature_residual(f::F, f_args::Fargs, f_hc, L, H, T) where {F,Fargs}
    
Helper function for updating θw, C, and the residual.
"""
@inline function temperature_residual(f::F, f_args::Fargs, f_hc, L, H, T) where {F,Fargs}
    θw = f(T, f_args...)
    C = f_hc(θw)
    Tres = T - (H - θw*L) / C
    return Tres, θw, C
end
"""
    DallAmico{Tprop,Tg,Tswrc<:VanGenuchten} <: SFCCFunction

Dall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.
"""
Base.@kwdef struct DallAmico{Tprop,Tg,Tswrc<:VanGenuchten} <: SFCCFunction
    prop::Tprop = SoilFreezingProperties()
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
@inline function (f::DallAmico)(T, θtot=f.prop.θtot, θsat=f.prop.θsat, θres=f.prop.θres, Tₘ=f.prop.Tₘ, α=f.swrc.α, n=f.swrc.n)
    let θsat = max(θtot, θsat),
        g = f.g,
        Lsl = f.prop.Lsl,
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
    DallAmicoSalt{Tprop,Tsc,TR,Tg,Tswrc} <: SFCCFunction

Freeze curve from Dall'Amico (2011) with saline freezing point depression.

Angelopoulos M, Westermann S, Overduin P, Faguet A, Olenchenko V, Grosse G, Grigoriev MN. Heat and salt flow in subsea permafrost
    modeled with CryoGRID2. Journal of Geophysical Research: Earth Surface. 2019 Apr;124(4):920-37.
"""
Base.@kwdef struct DallAmicoSalt{Tprop,Tsc,TR,Tg,Tswrc} <: SFCCFunction
    prop::Tprop = SoilFreezingProperties()
    saltconc::Tsc = 890.0u"mol/m^3" # salt concentration
    R::TR = 8.314459u"J/K/mol" #[J/K mol] universal gas constant
    g::Tg = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tswrc = VanGenuchten() # soil water retention curve
end
# DallAmico freeze curve with salt
function (f::DallAmicoSalt)(T,θsat,θtot,L,θres=f.θres,Tₘ=f.Tₘ,saltconc=f.saltconc,α=f.swrc.α,n=f.swrc.n)
    let θsat = max(θtot, θsat),
        g = f.prop.g,
        R = f.R,
        ρw = f.prop.ρw,
        m = 1-1/n,
        Tₘ = normalize_temperature(Tₘ),
        T = normalize_temperature(T),
        # freezing point depression based on salt concentration
        Tstar = Tₘ + Tₘ/L*(-R * saltconc * Tₘ),
        ψ₀ = f.swrc(inv, θtot, θsat, θres, α, n),
        # water matric potential
        ψ = ψ₀ + L / (ρw * g * Tstar) * (T - Tstar) * heaviside(Tstar-T),
        ψ = IfElse.ifelse(ψ < 0.0, ψ, 0.0);
        # van Genuchten evaulation to get θw
        return f.swrc(ψ,θres,θsat,α,n)
    end
end
"""
    McKenzie <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
Base.@kwdef struct McKenzie{Tprop,Γ} <: SFCCFunction
    prop::Tprop = SoilFreezingProperties()
    γ::Γ = 0.1u"K"
end
function (f::McKenzie)(T, θtot=f.prop.θtot, θsat=f.prop.θsat, θres=f.prop.θres, Tₘ=f.prop.Tₘ, γ=f.γ)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θsat = max(θtot, θsat);
        return IfElse.ifelse(T <= Tₘ, θres + (θsat-θres)*exp(-((T-Tₘ)/γ)^2), θtot)
    end
end
"""
    Westermann <: SFCCFunction

Westermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of
    wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,
    https://doi.org/10.5194/tc-5-945-2011, 2011. 
"""
Base.@kwdef struct Westermann{Tprop,Δ} <: SFCCFunction
    prop::Tprop = SoilFreezingProperties()
    δ::Δ = 0.1u"K"
end
function (f::Westermann)(T, θtot=f.prop.θtot, θsat=f.prop.θsat, θres=f.prop.θres, Tₘ=f.prop.Tₘ, δ=f.δ)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θsat = max(θtot, θsat);
        return IfElse.ifelse(T <= Tₘ, θres - (θsat-θres)*(δ/(T-Tₘ-δ)), θtot)
    end
end
