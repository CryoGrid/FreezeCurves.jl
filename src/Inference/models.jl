"""
    swvmodel(swp::SoilWaterVolume; sat_dispersion=10.0, por_dispersion=10.0, res_dispersion=10.0)

Probabilistic model of the porous water volume. Considers saturation `sat`, porosity `por`, and the residual fraction `res`
as random variables, returning `sat`, `θsat = por`, and `θres = res*por` as generated quantities.
"""
@model function swvmodel(vol::SoilWaterVolume; sat_dispersion=10.0, por_dispersion=10.0, res_dispersion=10.0)
    sat_mean ~ betaprior(1.0, sat_dispersion)
    sat_disp ~ Exponential(sat_dispersion-1)
    por_mean ~ betaprior(vol.θsat, por_dispersion)
    por_disp ~ Exponential(por_dispersion-1)
    sat ~ betaprior(sat_mean, 1 + sat_disp)
    por ~ betaprior(por_mean, 1 + por_disp)
    res ~ betaprior(vol.θres / vol.θsat, res_dispersion)
    θres = res*por
    θsat = por
    return (; sat, θsat, θres)
end

"""
    sfccpriors(model::SFCCModel)

Defines default hyperpriors and/or hyperparameters for `model`. Implementations should return a `NamedTuple` which
can then be passed to `sfccmodel`.
"""
sfccpriors(model::SFCCModel) = error("missing implementation for $(typeof(model))")
"""
    sfccmodel(model::SFCCModel{Tcfg,Tfc}, T::AbstractVector, priors, args...; kwargs...) where {Tcfg,Tfc}

Create a `Turing` probabilistic model for the given soil freezing characteristic curve and model configuration.
All implementations should accept temperature `T` (in °C) as an input argument and define a predicted variable
`θ` which represents the liquid water content.
"""
sfccmodel(model::SFCCModel, T::AbstractVector, priors, args...; kwargs...) = error("missing implementation for $(typeof(model))")

# Simple, two stage models

sfccpriors(::SFCCModel{<:Westermann,IsoNormal}) = (
    logδ = Normal(0,2),
    Tₘ = truncated(Normal(0,0.5), -Inf, 0),
    σ = Exponential(0.5),
    res_dispersion = 10.0,
    por_dispersion = 10.0,
    sat_dispersion = 10.0,
)
sfccpriors(::SFCCModel{<:McKenzie,IsoNormal}) = (
    logγ = Normal(0,2),
    Tₘ = truncated(Normal(0,0.5), -Inf, 0),
    σ = Exponential(0.5),
    res_dispersion = 10.0,
    por_dispersion = 10.0,
    sat_dispersion = 10.0,
)
const param_name_map = Dict(
    Westermann => :δ,
    McKenzie => :γ,
)
@model function sfccmodel(
    model::SFCCModel{Tfc,IsoNormal},
    T::AbstractVector,
    priors=sfccpriors(model),
) where {Tfc<:Union{McKenzie,Westermann}}
    fc = model.sfcc.f
    pname = param_name_map[Tfc]
    logp ~ NamedDist(dist, pname)
	Tₘ ~ priors.Tₘ
    p = exp(logp)
    sat_dispersion, por_dispersion, res_dispersion = priors.sat_dispersion, priors.por_dispersion, priors.res_dispersion
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc); sat_dispersion, por_dispersion, res_dispersion)
    θ_pred = fc.(T, sat; θsat, θres, Tₘ, pname => p)
    # truncated normal likelihood
    σ ~ priors.σ
    θ ~ arraydist(truncated.(Normal.(θ_pred,σ), 0, 1))
    return θ_pred
end

sfccpriors(::SFCCModel{<:DallAmico,IsoNormal}) = (
    logα = Normal(0,2),
    logn = Normal(0,2),
    Tₘ = truncated(Normal(0,0.5), -Inf, 0),
    σ = Exponential(0.5),
    res_dispersion = 10.0,
    por_dispersion = 10.0,
    sat_dispersion = 10.0,
)
@model function sfccmodel(
    model::SFCCModel{<:DallAmico,IsoNormal},
    T::AbstractVector,
    priors=sfccpriors(model),
)
    fc = model.sfcc.f
	logα ~ priors.logα
	logn ~ priors.logn
	Tₘ ~ priors.Tₘ
    α = exp(logα)
    n = 1 + exp(logn)
    sat_dispersion, por_dispersion, res_dispersion = priors.sat_dispersion, priors.por_dispersion, priors.res_dispersion
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc); sat_dispersion, por_dispersion, res_dispersion)
    θ_pred = fc.(T, sat; θsat, θres, Tₘ, α, n)
    # truncated normal likelihood
    σ ~ priors.σ
    θ ~ arraydist(truncated.(Normal.(θ_pred,σ), 0, 1))
    return θ_pred
end

sfccpriors(::SFCCModel{<:PainterKarra,IsoNormal}) = (
    logα = Normal(0,2),
    logn = Normal(0,2),
    Tₘ = truncated(Normal(0,0.5), -Inf, 0),
    σ = Exponential(0.5),
    β = Exponential(1.0),
    ω₀ = Beta(1,1),
    res_dispersion = 10.0,
    por_dispersion = 10.0,
    sat_dispersion = 10.0,
)
@model function sfccmodel(
    model::SFCCModel{<:PainterKarra,IsoNormal},
    T::AbstractVector,
    priors=sfccpriors(model),
)
    fc = model.sfcc.f
	logα ~ priors.logα
	logn ~ priors.logn
	Tₘ ~ priors.Tₘ
    β ~ priors.β
    ω₀ ~ priors.ω₀
    ω = ω₀/β
    α = exp(logα)
    n = 1 + exp(logn)
    sat_dispersion, por_dispersion, res_dispersion = priors.sat_dispersion, priors.por_dispersion, priors.res_dispersion
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc); sat_dispersion, por_dispersion, res_dispersion)
    θ_pred = fc.(T, sat; θsat, θres, Tₘ, β, ω, α, n)
    # truncated normal likelihood
    σ ~ priors.σ
    θ ~ arraydist(truncated.(Normal.(θ_pred,σ), 0, 1))
    return θ_pred
end
