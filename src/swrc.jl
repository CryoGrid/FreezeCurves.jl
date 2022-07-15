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
Base.inv(f::SWRCFunction) = (args...) -> f(inv, args...)
"""
    VanGenuchten <: SWRCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.
"""
Base.@kwdef struct VanGenuchten{Tα,Tn} <: SWRCFunction
    α::Tα = 1.0u"1/m"
    n::Tn = 2.0
end
function (f::VanGenuchten)(ψ, θsat, θres, α=f.α, n=f.n)
    let m = 1-1/n;
        IfElse.ifelse(ψ <= zero(ψ), θres + (θsat - θres)*(1 + (-α*ψ)^n)^(-m), θsat)
    end
end
function (f::VanGenuchten)(::typeof(inv), θ, θsat, θres, α=f.α, n=f.n)
    let m = 1-1/n;
        IfElse.ifelse(θ < θsat, -1/α*(((θ-θres)/(θsat-θres))^(-1/m)-1.0)^(1/n), zero(1/α))
    end
end
