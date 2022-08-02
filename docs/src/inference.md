## Estimating freeze curve parameters

Most soil freezing characteristic curves have one or more parameters which determine the shape of the curve. These parameters are typically determined empirically and can vary considerably between different types of soil.

FreezeCurves.jl facilitates inference of these parameters via the [Turing.jl](https://turing.ml) probabilistic programming language. In order to keep FreezeCurves.jl as lightweight as possible, `Turing` is not included as a dependncy by default. It must be installed separately (`Pkg.add("Turing")`) to expose the `FreezeCurves.Inference` module.

We can demonstrate this with an idealized test case where simply corrupt the "true" liquid water content values with isotropic Gaussian noise. Here we fit the van Genuchten parameters for the Dall'Amico (2011) freezing characteristic curve using the Bayesian probabilistic freeze curve model pre-sepcified in the `Inference` module. This model assumes that the measurement errors are normally distributed, which matches our idealized test case.

```julia
using Turing
using FreezeCurves
using FreezeCurves.Inference

import Random

rng = Random.MersenneTwister(1234)

fc = DallAmico(swrc=VanGenuchten(α=0.1u"1/m", n=1.8))
Trange = vcat(-5.0u"°C":0.1u"K":-0.11u"°C", -0.1u"°C":0.001u"K":0.0u"°C")
θtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))
sfcc_model = SFCCModel(fc) # from FreezeCurves.Inference
m = sfcc_model(Trange, θtrue) # condition on data
# draw 1,000 samples using the No U-Turn Sampler; gradients are computed automatically by Turing using forward-mode automatic differentiation (ForwardDiff.jl).
chain = sample(rng, m, NUTS(), 1000)
display(chain)
# output
Chains MCMC chain (1000×18×1 Array{Float64, 3}):

Iterations        = 501:1:1500
Number of chains  = 1
Samples per chain = 1000
Wall duration     = 5.07 seconds
Compute duration  = 5.07 seconds
parameters        = logα, logn, Tₘ, θfrac, θres, σ
internals         = lp, n_steps, is_accept, acceptance_rate, log_density, hamiltonian_energy, hamiltonian_energy_error, max_hamiltonian_energy_error, tree_depth, numerical_error, step_size, nom_step_size

Summary Statistics
  parameters      mean       std   naive_se      mcse        ess      rhat   ess_per_sec 
      Symbol   Float64   Float64    Float64   Float64    Float64   Float64       Float64 

        logα   -2.2728    0.0630     0.0020    0.0060    67.7835    1.0039       13.3590
        logn   -0.2151    0.0513     0.0016    0.0039   145.5487    1.0033       28.6852
          Tₘ   -0.0051    0.0036     0.0001    0.0004    45.3223    1.0085        8.9323
       θfrac    0.4869    0.0049     0.0002    0.0004    80.6756    1.0039       15.8998
        θres    0.0040    0.0032     0.0001    0.0002   326.9962    1.0007       64.4455
           σ    0.0187    0.0012     0.0000    0.0001   243.4807    1.0007       47.9859

Quantiles
  parameters      2.5%     25.0%     50.0%     75.0%     97.5% 
      Symbol   Float64   Float64   Float64   Float64   Float64 

        logα   -2.3786   -2.3181   -2.2779   -2.2364   -2.1379
        logn   -0.3026   -0.2513   -0.2180   -0.1816   -0.1120
          Tₘ   -0.0145   -0.0074   -0.0046   -0.0023   -0.0002
       θfrac    0.4778    0.4836    0.4870    0.4901    0.4962
        θres    0.0001    0.0014    0.0034    0.0060    0.0117
           σ    0.0167    0.0179    0.0186    0.0195    0.0214
```

This implementation uses the log-transform of the van Genuchten parameters which, while not strictly necessary, has a few positive effects:

1) It transforms the shape parameters, which have bounded support, into unconstrained space allowing for simple, unconstrained Gaussian priors.

2) It linearizes the otherwise non-linear relationship between α and n.

3) It smooths the partial derivatives of the objective function w.r.t the parameters which can improve the performance of gradient-based sampling and optimization methods (NUTS is a variant of Hamiltonian Monte Carlo which uses the gradient).

We can manually obtain the mean estimates of original parameters as `α = mean(exp.(Array(group(chain, :logα)))) ≈ 0.103` and `n = mean(1 .+ exp.(Array(group(chain, :logn)))) ≈ 1.808`, which are quite close to the "true" values as expected.

Real measurement data can be used by simply replacing `Trange` and `θtrue` in this example with equal-length vectors of temperature and water content measurements. In cases where the isotropic Gaussian error assumption does not hold, the `SFCCModel` interface can be extended to use custom model implementations.

If you aren't interested in the full posterior distribution, Turing also provides a very convenient interface to get a *maximum a posteriori* (MAP) estimate (i.e. a mode of the posterior) using [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl):

```julia
using Optim

optimize(m, MAP(), LBFGS())
# output
ModeResult with maximized lp of 394.29
6-element Named Vector{Float64}
A      │ 
───────┼────────────
:logα  │    -2.28583
:logn  │    -0.24065
:Tₘ    │ -0.00340876
:θfrac │    0.490879
:θres  │ 2.05868e-21
:σ     │    0.018277
```

Alternatively, one can ignore the prior entirely and just get a *maximum likelihood estimate*:

```julia
# here we use the common LBFGS optimizer; see Optim.jl docs for more options
optimize(m, MLE(), LBFGS())
# output
ModeResult with maximized lp of 394.29
6-element Named Vector{Float64}
A      │ 
───────┼────────────
:logα  │    -2.28583
:logn  │    -0.24065
:Tₘ    │ -0.00340876
:θfrac │    0.490879
:θres  │ 8.58955e-23
:σ     │    0.018277
```
