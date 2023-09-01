"""
    IsoNormalVWC{T} <: SFCCLikelihood

Represents a (truncated) isotropic Gaussian distributed observation model for volumetric water content.
This seems to be generally preferabe to the Beta likelihood both numerically and in terms of interpretability.
"""
Base.@kwdef struct IsoNormalVWC{T} <: SFCCLikelihood
    σ_mean::T = 0.1
end
sfccpriors(m::IsoNormalVWC) = (
    σ = Exponential(m.σ_mean),
)
@model function sfcclikelihood(model::IsoNormalVWC, θ_pred, priors)
    # truncated normal likelihood
    σ ~ priors.σ
    θ ~ arraydist(truncated.(Normal.(θ_pred,σ), 0, 1))
    return θ
end

"""
    BetaVWC{T} <: SFCCLikelihood

Represents a Beta distribtued observation model for volumetric water content.
"""
Base.@kwdef struct BetaVWC{T} <: SFCCLikelihood
    dispersion_mean::T = 100.0
end
sfccpriors(m::BetaVWC) = (
    vwc_dispersion = truncated(Exponential(m.dispersion_mean), 1, Inf),
)
@model function sfcclikelihood(model::BetaVWC, θ_pred, priors)
    # truncated normal likelihood
    vwc_dispersion ~ priors.vwc_dispersion
    μ = θ_pred
    ϕ = vwc_dispersion
    # take max w/ sqrt(eps()) to prevent zero values
    a = max.(μ.*ϕ, sqrt(eps()))
    b = max.((1 .- μ).*ϕ, sqrt(eps()))
    θ ~ arraydist(Beta.(a, b))
    return θ
end

struct ExactTemperatures <: TemperatureMeasurementModel end
sfccpriors(::ExactTemperatures) = (;)
@model temperature_measurement_model(::ExactTemperatures, T_obs::AbstractVector, priors) = T_obs

Base.@kwdef struct IsoNormalT <: TemperatureMeasurementModel
    obs_noise = 0.1
end
sfccpriors(m::IsoNormalT) = (
    σ = Exponential(m.obs_noise),
)
@model function temperature_measurement_model(m::IsoNormalT, T_obs::AbstractVector, priors)
    T_μ ~ MvNormal(T_obs, 10*m.obs_noise)
    σ ~ priors.σ
    T ~ MvNormal(T_μ, σ)
    return T
end

sfccpriors(vol::SoilWaterVolume; sat_mean=1.0, sat_dispersion=10.0, por_dispersion=10.0, res_dispersion=10.0) = (
    sat = betaprior(sat_mean, sat_dispersion),
    por = betaprior(vol.θsat, por_dispersion),
    res = betaprior(vol.θres / vol.θsat, res_dispersion),
)
"""
    swvmodel(swp::SoilWaterVolume, priors)

Probabilistic model of the porous water volume. Considers saturation `sat`, porosity `por`, and the residual fraction `res`
as random variables, returning `sat`, `θsat = por`, and `θres = res*por` as generated quantities.
"""
@model function swvmodel(::SoilWaterVolume, priors)
    sat ~ priors.sat
    por ~ priors.por
    res ~ priors.res
    θres = res*por
    θsat = por
    return (; sat, θsat, θres)
end

"""
    sfccmodel(model::SFCCModel{Tcfg,Tfc}, T::AbstractVector, priors, args...; kwargs...) where {Tcfg,Tfc}

Create a `Turing` probabilistic model for the given soil freezing characteristic curve and model configuration.
All implementations should accept temperature `T` (in °C) as an input argument and define a predicted variable
`θ` which represents the liquid water content.
"""
sfccmodel(model::SFCCModel, T::AbstractVector, priors, args...; kwargs...) = error("missing implementation for $(typeof(model))")

sfccpriors(m::SFCCModel{<:Westermann}) = (
    logδ = Normal(0,2),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
sfccpriors(m::SFCCModel{<:McKenzie}) = (
    logγ = Normal(0,2),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
sfccpriors(m::SFCCModel{<:Hu2020}) = (
    b = Normal(0,2),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{Tfc},
    T_obs::AbstractVector;
    priors=(;),
) where {Tfc<:Union{McKenzie,Westermann}}
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
    pname = if Tfc <: McKenzie
        :γ
    else
        :δ
    end
    logp ~ NamedDist(dist, pname)
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    p = exp(logp)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θ_pred = fc.(T, sat; θsat, θres, pname => p)
    @submodel sfcclikelihood(model.lik, θ_pred, priors.lik)
    return θ_pred
end

sfccpriors(m::SFCCModel{<:Hu2020}) = (
    b = Beta(1,2),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:Hu2020},
    T_obs::AbstractVector;
    priors=(;),
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
    b ~ priors.b
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θ_pred = fc.(T, sat; θsat, θres, b)
    @submodel sfcclikelihood(model.lik, θ_pred, priors.lik)
    return θ_pred
end

sfccpriors(m::SFCCModel{<:DallAmico}) = (
    logα = Normal(0,2),
    logn = Normal(0,2),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:DallAmico},
    T_obs::AbstractVector,
    priors=(;),
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
	logα ~ priors.logα
	logn ~ priors.logn
    α = exp(logα)
    n = 1 + exp(logn)
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θ_pred = fc.(T, sat; θsat, θres, α, n)
    @submodel sfcclikelihood(model.lik, θ_pred, priors.lik)
    return θ_pred
end

sfccpriors(m::SFCCModel{<:DallAmicoSalt}) = (
    logα = Normal(0,2),
    logn = Normal(0,2),
    saltconc = Exponential(200.0),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:DallAmicoSalt},
    T_obs::AbstractVector,
    priors=(;),
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
	logα ~ priors.logα
	logn ~ priors.logn
    saltconc ~ prior.saltconc
    α = exp(logα)
    n = 1 + exp(logn)
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θ_pred = fc.(T, sat; θsat, θres, α, n, saltconc)
    @submodel sfcclikelihood(model.lik, θ_pred, priors.lik)
    return θ_pred
end

sfccpriors(m::SFCCModel{<:PainterKarra}) = (
    logα = Normal(0,2),
    logn = Normal(0,2),
    β = Exponential(1.0),
    ω₀ = Beta(1,1),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:PainterKarra},
    T_obs::AbstractVector,
    priors=(;),
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
	logα ~ priors.logα
	logn ~ priors.logn
    β ~ priors.β
    ω₀ ~ priors.ω₀
    ω = ω₀/β
    α = exp(logα)
    n = 1 + exp(logn)
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θ_pred = fc.(T, sat; θsat, θres, β, ω, α, n)
    @submodel sfcclikelihood(model.lik, θ_pred, priors.lik)
    return θ_pred
end
