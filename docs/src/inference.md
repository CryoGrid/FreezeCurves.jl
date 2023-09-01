## Estimating freeze curve parameters

Most soil freezing characteristic curves have one or more parameters which determine the shape of the curve. These parameters are typically determined empirically and can vary considerably between different types of soil.

FreezeCurves.jl facilitates inference of these parameters via the [Turing.jl](https://turing.ml) probabilistic programming language. In order to keep FreezeCurves.jl as lightweight as possible, `Turing` is not included as a dependncy by default. It must be installed separately (`Pkg.add("Turing")`) for this functionality to be enabled.

We can demonstrate this with an idealized test case where simply corrupt the "true" liquid water content values with isotropic Gaussian noise. Here we fit the van Genuchten parameters for the Dall'Amico (2011) freezing characteristic curve using the Bayesian probabilistic freeze curve model pre-sepcified in `FreezeCurves`. This model assumes that the measurement errors are normally distributed, which matches our idealized test case.

```julia
using Turing
using FreezeCurves

import Random

rng = Random.MersenneTwister(1234)

fc = DallAmico(swrc=VanGenuchten(α=0.1u"1/m", n=1.8))
Trange = vcat(-5.0u"°C":0.1u"K":-0.11u"°C", -0.1u"°C":0.001u"K":0.0u"°C")
θtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))
sfcc_model = SFCCModel(fc)
m = sfcc_model(Trange, θtrue) # condition on data
# draw 1,000 samples using the No U-Turn Sampler w/ 500 adaptation steps and 85% target acceptance rate; gradients are computed automatically by Turing using forward-mode automatic differentiation (ForwardDiff.jl).
chain = sample(rng, m, NUTS(500,0.85), 1_000)
display(chain)
```
Output:
```
Chains MCMC chain (1000×19×1 Array{Float64, 3}):

Iterations        = 501:1:1500
Number of chains  = 1
Samples per chain = 1000
Wall duration     = 10.09 seconds
Compute duration  = 10.09 seconds
parameters        = logα, logn, Tₘ, sat, por, res, σ
internals         = lp, n_steps, is_accept, acceptance_rate, log_density, hamiltonian_energy, hamiltonian_energy_error, max_hamiltonian_energy_error, tree_depth, numerical_error, step_size, nom_step_size

Summary Statistics
  parameters      mean       std   naive_se      mcse        ess      rhat   ess_per_sec 
      Symbol   Float64   Float64    Float64   Float64    Float64   Float64       Float64 

        logα   -2.2574    0.0575     0.0018    0.0033   381.2072    1.0010       37.7732
        logn   -0.2221    0.0494     0.0016    0.0026   456.1224    0.9990       45.1964
          Tₘ   -0.0034    0.0027     0.0001    0.0001   371.9088    0.9994       36.8518
         sat    0.9889    0.0101     0.0003    0.0004   437.2432    1.0006       43.3257
         por    0.4949    0.0052     0.0002    0.0002   368.9660    1.0042       36.5602
         res    0.0074    0.0063     0.0002    0.0002   631.1171    1.0044       62.5364
           σ    0.0193    0.0012     0.0000    0.0000   787.7707    0.9992       78.0589

Quantiles
  parameters      2.5%     25.0%     50.0%     75.0%     97.5% 
      Symbol   Float64   Float64   Float64   Float64   Float64 

        logα   -2.3575   -2.2960   -2.2622   -2.2221   -2.1284
        logn   -0.3147   -0.2585   -0.2219   -0.1887   -0.1282
          Tₘ   -0.0100   -0.0050   -0.0028   -0.0012   -0.0001
         sat    0.9621    0.9848    0.9919    0.9962    0.9997
         por    0.4851    0.4915    0.4947    0.4980    0.5053
         res    0.0002    0.0025    0.0058    0.0107    0.0238
           σ    0.0171    0.0184    0.0192    0.0200    0.0218
```

This implementation uses the log-transform of the van Genuchten parameters which, while not strictly necessary, has a few positive effects:

1) It transforms the shape parameters, which have bounded support, into unconstrained space allowing for simple, unconstrained Gaussian priors.

2) It linearizes the otherwise non-linear relationship between α and n.

3) It smooths the partial derivatives of the objective function w.r.t the parameters which can improve the performance of gradient-based sampling and optimization methods (NUTS is a variant of Hamiltonian Monte Carlo which uses the gradient).

We can manually obtain the mean estimates of original parameters as `α = mean(exp.(Array(group(chain, :logα)))) ≈ 0.105` and `n = mean(1 .+ exp.(Array(group(chain, :logn)))) ≈ 1.802`, which are quite close to the "true" values as expected.

Real measurement data can be used by simply replacing `Trange` and `θtrue` in this example with equal-length vectors of temperature and water content measurements. In cases where the isotropic Gaussian error assumption does not hold, the `SFCCModel` interface can be extended to use custom model implementations.

If you aren't interested in the full posterior distribution, Turing also provides a very convenient interface to get a *maximum a posteriori* (MAP) estimate (i.e. a mode of the posterior) using [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl):

```julia
using Optim

optimize(m, MAP(), LBFGS())
```
Output:
```
ModeResult with maximized lp of 394.16
7-element Named Vector{Float64}
A     │ 
──────┼─────────────
:logα │      -2.3055
:logn │    -0.228226
:Tₘ   │ -1.45969e-52
:sat  │          1.0
:por  │     0.496235
:res  │  4.83361e-31
:σ    │     0.018805
```

Alternatively, one can ignore the prior entirely and just get a *maximum likelihood estimate*:

```julia
# here we use the common LBFGS optimizer; see Optim.jl docs for more options
optimize(m, MLE(), LBFGS())
```
Output:
```
ModeResult with maximized lp of 394.16
7-element Named Vector{Float64}
A     │ 
──────┼─────────────
:logα │      -2.3055
:logn │    -0.228226
:Tₘ   │ -6.85422e-72
:sat  │          1.0
:por  │     0.496235
:res  │  4.64697e-43
:σ    │     0.018805
```

Note how the MLE vs MAP results are almost identical in this case since we have a (relatively) large generated dataset and the default priors are only very weakly informative.
