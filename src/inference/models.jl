"""
    IsoNormalVWC{T} <: SFCCLikelihood

Represents a (truncated) isotropic Gaussian distributed observation model for volumetric water content.
This seems to be generally preferabe to the Beta likelihood both numerically and in terms of interpretability.
"""
Base.@kwdef struct IsoNormalVWC{T} <: SFCCLikelihood
    σ_mean::T = 0.01
end
sfccpriors(m::IsoNormalVWC) = (
    σ = Exponential(m.σ_mean),
)
@model function sfcclikelihood(model::IsoNormalVWC, pred, priors)
    # truncated normal likelihood
    σ ~ priors.σ
    θ ~ arraydist(truncated.(Normal.(pred.θw, σ), 0, 1))
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
@model function sfcclikelihood(model::BetaVWC, pred, priors)
    # truncated normal likelihood
    vwc_dispersion ~ priors.vwc_dispersion
    μ = pred.θw
    ϕ = vwc_dispersion
    # take max w/ sqrt(eps()) to prevent zero values
    a = max.(μ.*ϕ, sqrt(eps()))
    b = max.((1 .- μ).*ϕ, sqrt(eps()))
    θ ~ arraydist(Beta.(a, b))
    return θ
end

"""
    BetaNormalVWC{T} <: SFCCLikelihood

Represents a Beta distribtued observation model for volumetric water content with
normally distribtued noise above the melting point.
"""
Base.@kwdef struct BetaNormalVWC{Ts1,Td1,Td2} <: SFCCLikelihood
    σ_sat_mean::Ts1 = 0.05
    dispersion_mean::Td1 = 99.0
    dispersion_scale::Td2 = 50.0
end
sfccpriors(m::BetaNormalVWC) = (
    σ_sat = Exponential(m.σ_sat_mean),
    vwc_dispersion = from_moments(LogNormal, m.dispersion_mean, m.dispersion_scale),
)
@model function sfcclikelihood(model::BetaNormalVWC, pred, priors)
    function vwcdist(θw, T, Tₘ, disp, σ_sat)
        if T >= Tₘ
            return Normal(θw, σ_sat)
        else
            μ = θw
            ϕ = 1 + disp
            # take max w/ sqrt(eps()) to prevent zero values
            a = max(μ*ϕ, sqrt(eps()))
            b = max((1 - μ)*ϕ, sqrt(eps()))
            return Beta(a, b)
        end
    end
    σ_sat ~ priors.σ_sat
    vwc_dispersion ~ priors.vwc_dispersion
    θ ~ arraydist(vwcdist.(pred.θw, pred.T, pred.Tₘ, vwc_dispersion, σ_sat))
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
@model function swvmodel(vol::SoilWaterVolume, priors)
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
    Tₘ = truncated(Normal(0,0.25), -Inf, 0),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
sfccpriors(m::SFCCModel{<:McKenzie}) = (
    logγ = Normal(0,2),
    Tₘ = truncated(Normal(0,0.25), -Inf, 0),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
sfccpriors(m::SFCCModel{<:Hu2020}) = (
    b = Normal(0,2),
    Tₘ = truncated(Normal(0,0.25), -Inf, 0),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{Tfc},
    T_obs::AbstractVector;
    priors=(;),
    temperature_measurement_model=temperature_measurement_model,
    swvmodel=swvmodel,
    sfcclikelihood=sfcclikelihood,
) where {Tfc<:Union{McKenzie,Westermann}}
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
    pname = if Tfc <: McKenzie
        :γ
    else
        :δ
    end
    logp ~ NamedDist(dist, pname)
	Tₘ ~ priors.Tₘ
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    p = exp(logp)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θw = fc.(T, sat; θsat, θres, Tₘ, pname => p)
    pred = (; θw, θsat, θres, T, Tₘ)
    @submodel sfcclikelihood(model.lik, pred, priors.lik)
    return pred
end

sfccpriors(m::SFCCModel{<:Hu2020}) = (
    b = Beta(1,2),
    Tₘ = truncated(Normal(0,0.25), -Inf, 0),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:Hu2020},
    T_obs::AbstractVector;
    priors=(;),
    temperature_measurement_model=temperature_measurement_model,
    swvmodel=swvmodel,
    sfcclikelihood=sfcclikelihood,
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
    b ~ priors.b
	Tₘ ~ priors.Tₘ
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    p = exp(logp)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θw = fc.(T, sat; θsat, θres, Tₘ, b)
    pred = (; θw, θsat, θres, T, Tₘ)
    @submodel sfcclikelihood(model.lik, pred, priors.lik)
    return pred
end

sfccpriors(m::SFCCModel{<:DallAmico}) = (
    logα = Normal(0,5),
    logn = Normal(0,2),
    Tₘ = truncated(Normal(0,0.25), -Inf, 0),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:DallAmico},
    T_obs::AbstractVector,
    priors=(;),
    temperature_measurement_model=temperature_measurement_model,
    swvmodel=swvmodel,
    sfcclikelihood=sfcclikelihood,
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
	logα ~ priors.logα
	logn ~ priors.logn
	Tₘ ~ priors.Tₘ
    α = exp(logα)
    n = 1 + exp(logn)
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θw = fc.(T, sat; θsat, θres, Tₘ, α, n)
    pred = (; θw, θsat, θres, T, Tₘ)
    @submodel sfcclikelihood(model.lik, pred, priors.lik)
    return pred
end

sfccpriors(m::SFCCModel{<:PainterKarra}) = (
    logα = Normal(0,5),
    logn = Normal(0,2),
    Tₘ = truncated(Normal(0,0.25), -Inf, 0),
    # β = Exponential(1.0),
    ω₀ = Beta(1,1),
    lik = sfccpriors(m.lik),
    meas = sfccpriors(m.meas),
    vol = sfccpriors(SoilWaterVolume(m.sfcc)),
)
@model function sfccmodel(
    model::SFCCModel{<:PainterKarra},
    T_obs::AbstractVector,
    priors=(;),
    temperature_measurement_model=temperature_measurement_model,
    swvmodel=swvmodel,
    sfcclikelihood=sfcclikelihood,
)
    priors = merge(sfccpriors(model), priors)
    fc = model.sfcc
	logα ~ priors.logα
	logn ~ priors.logn
	Tₘ ~ priors.Tₘ
    # β ~ priors.β
    β = fc.β
    ω₀ ~ priors.ω₀
    ω = ω₀/β
    α = exp(logα)
    n = 1 + exp(logn)
    T = @submodel temperature_measurement_model(model.meas, T_obs, priors.meas)
    sat, θsat, θres = @submodel swvmodel(SoilWaterVolume(fc), priors.vol)
    θw = fc.(T, sat; θsat, θres, Tₘ, β, ω, α, n)
    pred = (; θw, θsat, θres, T, Tₘ)
    @submodel sfcclikelihood(model.lik, pred, priors.lik)
    return pred
end

function from_moments(::Type{LogNormal}, mean, stddev)
    var = stddev^2
    μ = log(mean / sqrt(var / mean^2 + 1))
    σ = sqrt(log(var / mean^2 + 1))
    return LogNormal(μ, σ)
end
