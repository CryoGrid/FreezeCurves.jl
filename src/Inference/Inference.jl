module Inference

using ..FreezeCurves
using ..FreezeCurves.Solvers

using ..Turing
using ..Turing.Distributions

export SFCCModel, sfccpriors

betaprior(mean, dispersion) = Beta(max(mean*dispersion,1), max((1-mean)*dispersion,1))

abstract type SFCCLikelihood end

abstract type TemperatureMeasurementModel end

"""
    SFCCModel{Tfc<:SFCCFunction,Tlik,Tsp,Tsolver<:SFCCSolver}

Represents a statistical model for the soil freezing characteristic curve specified by `Tfc` with configuration `Tcfg`.
"""
struct SFCCModel{Tfc<:SFCCFunction,Tlik<:SFCCLikelihood,Tmeas<:TemperatureMeasurementModel,Tsp,Tsolver<:SFCCSolver}
    sfcc::SFCC{Tfc,Tsolver}
    lik::Tlik
    meas::Tmeas
    sp::Tsp
    function SFCCModel(fc::SFCCFunction, sp=nothing; likelihood=IsoNormalVWC(), meas=ExactTemperatures(), solver::SFCCSolver=SFCCNewtonSolver())
        fc = ustrip(fc)
        return new{typeof(fc),typeof(likelihood),typeof(meas),typeof(sp),typeof(solver)}(SFCC(fc, solver), likelihood, meas, sp)
    end
end
"""
    (model::SFCCModel)(T::AbstractVector, θ::AbstractVector, args...; kwargs...)

Constructs a `Turing` model with `sfccmodel` given temperatures `T`, conditioned on `θ`.
The unconditional model (with `θ` sampled from the model) can be obtained using `Turing.decondition`.
"""
function (model::SFCCModel)(T::AbstractVector, θ::AbstractVector, args...; kwargs...)
    @assert length(T) == length(θ)
    T = ustrip.(u"°C", T)
    m = sfccmodel(model, T, args...; kwargs...)
    m_cond = if isa(model.meas, ExactTemperatures)
        Turing.condition(m; θ)
    else
        Turing.condition(m; θ, T)
    end
    return m_cond
end

"""
    preprocess_data(T, θ; T_min=-30.0, T_max=0.0)

Utility method for preprocessing (paired) temperature and volumetric water content measurements.
Drops all measurements with `missing` values, invalid water contents (`< 0` or `> 1`), or temperatures (°C)
outside of a given range.
"""
function preprocess_data(T, θ; T_min=-30.0, T_max=0.0)
    @assert length(T) == length(θ)
    is_invalid = ismissing.(T) .| ismissing.(θ) .| (θ .< 0.0) .| (θ  > 1.0) .| (T .< T_min) .| (T .> T_max)
    T = ifelse.(is_invalid, missing, T_all) |> skipmissing |> collect
    θ = ifelse.(is_invalid, missing, θ_all) |> skipmissing |> collect
    return (; T, θ)
end

"""
    sfccpriors(model; kwargs...)

Defines default priors and/or hyperparameters for `model` which may be a `SFCCModel` or one of its
components (e.g. `SFCCLikelihood`). Implementations should return a `NamedTuple` where the names
match the model or submodel variables.
"""
function sfccpriors end

include("models.jl")

end
