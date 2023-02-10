var documenterSearchIndex = {"docs":
[{"location":"api/Solvers/#Solvers","page":"Solvers","title":"Solvers","text":"","category":"section"},{"location":"api/Solvers/","page":"Solvers","title":"Solvers","text":"The Solvers sub-module provides an interface for performing non-linear root-finding, e.g. solving implicit freeze curve formulations or recovering temperature given enthalpy, on the soil freeze characteristic curves implemented in this package. Note that this is not to be confused with parameter inference which is handled by a separate sub-module.","category":"page"},{"location":"api/Solvers/","page":"Solvers","title":"Solvers","text":"DocTestSetup = quote\n    using FreezeCurves\n    using FreezeCurves.Solvers\nend","category":"page"},{"location":"api/Solvers/","page":"Solvers","title":"Solvers","text":"Modules = [FreezeCurves.Solvers]\nPrivate = false\nOrder = [:type, :function, :macro]","category":"page"},{"location":"api/Solvers/#FreezeCurves.Solvers.SFCCInverseEnthalpyObjective","page":"Solvers","title":"FreezeCurves.Solvers.SFCCInverseEnthalpyObjective","text":"SFCCInverseEnthalpyObjective{TF,Tkwargs<:NamedTuple,Thc,TL,TH,Tsat} <: AbstractSFCCObjective\n\nOptimization objective for finding a temperature, T, wich resolves the conservation equation:\n\n0 = (H - Lθ)C - T\n\ngiven a fixed enthalpy value H and freeze curve arguments f_args.\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.SFCCNewtonSolver","page":"Solvers","title":"FreezeCurves.Solvers.SFCCNewtonSolver","text":"SFCCNewtonSolver <: SFCCSolver\n\nFast, specialized implementation of Newton's method with backtracking line search for resolving the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid jumping over the solution. This prevents convergence issues that arise due to discontinuities and strong non-linearity in most common soil freeze curves.\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.SFCCPreSolver","page":"Solvers","title":"FreezeCurves.Solvers.SFCCPreSolver","text":"SFCCPreSolver{TCache} <: SFCCSolver\n\nA fast SFCC \"solver\" which pre-builds an interpolant for the freeze curve in terms of enthalpy, θ(H). Note that this solver is only valid when all freeze curve parameters are held constant and will produce incorrect results otherwise.\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.heatcapacity-Tuple{Number, Number}","page":"Solvers","title":"FreezeCurves.Solvers.heatcapacity","text":"heatcapacity(c_frozen::Number, c_thawed::Number; θtot=1.0)\n\nPiecewise constant heat capacity function:\n\nf(θw) =\n    begincases\n    c_f  θw  θtot \n    c_t  textotherwise\n    endcases\n\n\n\n\n\n","category":"method"},{"location":"api/Solvers/#FreezeCurves.Solvers.heatcapacity-Tuple{Number}","page":"Solvers","title":"FreezeCurves.Solvers.heatcapacity","text":"heatcapacity(c::Number)\n\nHeat capacity function f(H,T,θw) = c where c is a pre-specified constant.\n\n\n\n\n\n","category":"method"},{"location":"api/Solvers/#FreezeCurves.Solvers.sfccsolve-Union{Tuple{return_all}, Tuple{FreezeCurves.Solvers.AbstractSFCCObjective, SFCCSolver, Any}, Tuple{FreezeCurves.Solvers.AbstractSFCCObjective, SFCCSolver, Any, Val{return_all}}} where return_all","page":"Solvers","title":"FreezeCurves.Solvers.sfccsolve","text":"sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀, ::Val{return_all}=Val{true}()) where {return_all}\n\nSolve the given objective obj using solver and initial guess x₀. If return_all=true, then sfccsolve should return a named tuple with at least the temperature solution T, the liquid water content θw, the heat capacity C, and the liquid water content deriviative ∂θw∂T defined. Solver-specific additional values may also be included.\n\n\n\n\n\n","category":"method"},{"location":"api/Solvers/#FreezeCurves.Solvers.sfccsolve-Union{Tuple{return_all}, Tuple{FreezeCurves.Solvers.SFCCInverseEnthalpyObjective, FreezeCurves.Solvers.SFCCNewtonSolver, Number}, Tuple{FreezeCurves.Solvers.SFCCInverseEnthalpyObjective, FreezeCurves.Solvers.SFCCNewtonSolver, Number, Val{return_all}}} where return_all","page":"Solvers","title":"FreezeCurves.Solvers.sfccsolve","text":"sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNewtonSolver, T₀::Number, ::Val{return_all}=Val{true}(); error_when_not_converged=true)\n\nSolves obj using the specialized Newton solver and returns the result. If return_all is true, a named tuple (;T, Tres, θw, C, itercount) is returned; otherwise (by default), only the temperature solution is returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#Index-of-public-API","page":"Index of public API","title":"Index of public API","text":"","category":"section"},{"location":"api/","page":"Index of public API","title":"Index of public API","text":"","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"api/FreezeCurves/#FreezeCurves","page":"FreezeCurves","title":"FreezeCurves","text":"","category":"section"},{"location":"api/FreezeCurves/","page":"FreezeCurves","title":"FreezeCurves","text":"Modules = [FreezeCurves]\nPrivate = false\nOrder = [:type, :function, :macro]","category":"page"},{"location":"api/FreezeCurves/#FreezeCurves.BrooksCorey","page":"FreezeCurves","title":"FreezeCurves.BrooksCorey","text":"BrooksCorey{Twp,Tψₛ,Tλ} <: SWRCFunction\n\nvan Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.     Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.DallAmico","page":"FreezeCurves","title":"FreezeCurves.DallAmico","text":"DallAmico{TFT,Tg,Tswrc<:SWRCFunction} <: SFCCFunction\n\nDall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.DallAmicoSalt","page":"FreezeCurves","title":"FreezeCurves.DallAmicoSalt","text":"DallAmicoSalt{TFT,Tsc,TR,Tg,Tswrc<:SWRCFunction} <: SFCCFunction\n\nFreeze curve from Dall'Amico (2011) with saline freezing point depression.\n\nAngelopoulos M, Westermann S, Overduin P, Faguet A, Olenchenko V, Grosse G, Grigoriev MN. Heat and salt flow in subsea permafrost     modeled with CryoGRID2. Journal of Geophysical Research: Earth Surface. 2019 Apr;124(4):920-37.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.Hu2020","page":"FreezeCurves","title":"FreezeCurves.Hu2020","text":"Hu2020{TFT,Tvol,Tb}  <: SFCCFunction\n\nSoil freezing characteristic curve formulation of Hu et al. 2020.\n\nHu G, Zhao L, Zhu X, Wu X, Wu T, Li R, Xie C, Hao J. Review of algorithms and parameterizations to determine unfrozen water content in frozen soil. Geoderma. 2020 Jun 1;368:114277.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.McKenzie","page":"FreezeCurves","title":"FreezeCurves.McKenzie","text":"McKenzie{TFT,Tvol,Tγ} <: SFCCFunction\n\nMcKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:     numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,     30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.PainterKarra","page":"FreezeCurves","title":"FreezeCurves.PainterKarra","text":"PainterKarra{TFT,Tβ,Tω,Tg,Tswrc<:SWRCFunction} <: SFCCFunction\n\nPainter SL, Karra S. Constitutive Model for Unfrozen Water Content in Subfreezing Unsaturated Soils. Vadose zone j. 2014 Apr;13(4):1-8.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.PowerLaw","page":"FreezeCurves","title":"FreezeCurves.PowerLaw","text":"PowerLaw{Tvol,Tα,Tβ} <: SFCCFunction\n\nCommonly used power law relation of Lovell (1957).\n\nLovell, C.: Temperature effects on phase composition and strength of partially frozen soil, Highway Research Board Bulletin, 168, 74–95, 1957.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SFCC","page":"FreezeCurves","title":"FreezeCurves.SFCC","text":"SFCC{F,S} <: FreezeCurve\n\nGeneric representation of the soil freezing characteristic curve along with a nonlinear solver for resolving the temperature-energy conservation law. The shape and parameters of the curve are determined by the implementation of SFCCFunction f.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SFCCFunction","page":"FreezeCurves","title":"FreezeCurves.SFCCFunction","text":"SFCCFunction\n\nBase type for a soil freeze characteristic curve (SFCC) function. Subtypes should be callable structs that implement the freeze curve and contain any necessary additional constants or configuration options. User-specified parameters can either be supplied in the struct or declared as model parameters via the variables method.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SFCCSolver","page":"FreezeCurves","title":"FreezeCurves.SFCCSolver","text":"SFCCSolver\n\nBase type representing non-linear solvers for implicit SFCC functions.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SWRCFunction","page":"FreezeCurves","title":"FreezeCurves.SWRCFunction","text":"SWRCFunction\n\nBase type for soil water retention curves (SWRC) which relate soil water matric potential to water content.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SoilFreezeThawProperties","page":"FreezeCurves","title":"FreezeCurves.SoilFreezeThawProperties","text":"SoilFreezeThawProperties{TTₘ,TLsl}\n\nStruct containing constants/parameters common to some or all SFCCs.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SoilFreezeThawProperties-Tuple{SFCCFunction}","page":"FreezeCurves","title":"FreezeCurves.SoilFreezeThawProperties","text":"SoilFreezeThawProperties(f::SFCCFunction)\n\nRetrieves the default SoilFreezeThawProperties from f; should be defined for all freeze curves.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.SoilWaterVolume","page":"FreezeCurves","title":"FreezeCurves.SoilWaterVolume","text":"SoilWaterVolume{Tρw,Tθres,Tθsat}\n\nStruct containing basic physical constants and propeties related to soil pore water.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SoilWaterVolume-Tuple{SFCCFunction}","page":"FreezeCurves","title":"FreezeCurves.SoilWaterVolume","text":"SoilWaterVolume(f::SFCCFunction)\n\nRetrieves the nested SoilWaterVolume from the SoilFreezeThawProperties of the freeze curve f.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.SoilWaterVolume-Tuple{SWRCFunction}","page":"FreezeCurves","title":"FreezeCurves.SoilWaterVolume","text":"SoilWaterVolume(f::SWRCFunction)\n\nGet the SoilWaterVolume defined for the given SWRCFunction f. Must be implemented for all SWRCFunction types. Default implementation retrieves the field water.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.VanGenuchten","page":"FreezeCurves","title":"FreezeCurves.VanGenuchten","text":"VanGenuchten{Tvol,Tα,Tn} <: SWRCFunction\n\nvan Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.     Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.Westermann","page":"FreezeCurves","title":"FreezeCurves.Westermann","text":"Westermann{TFT,Tvol,Tδ} <: SFCCFunction\n\nWestermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of     wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,     https://doi.org/10.5194/tc-5-945-2011, 2011. \n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.inflectionpoint-Tuple{SFCCFunction}","page":"FreezeCurves","title":"FreezeCurves.inflectionpoint","text":"inflectionpoint(f::SFCCFunction)\n\nReturns the analytical inflection point (i.e. where ∂²θ/∂T^2 = 0), if available.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.inflectionpoint-Tuple{SWRCFunction}","page":"FreezeCurves","title":"FreezeCurves.inflectionpoint","text":"inflectionpoint(f::SWRCFunction)\n\nReturns the analytical solution for the inflection point (i.e. where ∂²θ/∂ψ² = 0) of the SWRC, if available.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#Unitful.ustrip-Tuple{Union{SFCCFunction, SWRCFunction}}","page":"FreezeCurves","title":"Unitful.ustrip","text":"ustrip(x)\n\nReconstructs the type or function x with all numerical quantities stripped of units.\n\n\n\n\n\n","category":"method"},{"location":"#FreezeCurves.jl","page":"Home","title":"FreezeCurves.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FreezeCurves is a lightweight Julia package to facilitate the study and analysis of soil freezing characteristics (Koopmans and Miller, 1966). The relationship between temperature and unfrozen water content in porous media (such as soils) is often highly nonlinear and plays a significant role in the analysis of freeze-thaw dynamics in science and engineering. One common application is in the geophysical modeling of permafrost, where having faithful representations of freeze-thaw processes is often paramount to accurately resolving long-term changes in the subsurface thermal regime.","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"index.md\",\"inference.md\",\"api/index.md\"]","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The soil freezing characteristic curve (SFCC) is typically a monotonic function which maps (usually subzero) temperatures [°C] to volumetric unfrozen/liquid water contents [m³/m³]. The SFCC is closely related to the soil-water retention curve (SWRC) which governs the similar non-linear realtionship between soil-water matric potential [m] saturation [m³/m³]. SWRCs are widely used in hydrological modeling, with the two most common formulations being the van Genuchten (1980) and Brooks-Corey (Brooks, 1965) models. SFCCs are less widely known, but numerous models have been proposed over the years (Kurylyk and Watanbe, 2013).","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package is intended to act as a living repository of SFCC/SWRC implementations which can then be fitted to data or consumed downstream by thermal or hydrological models such as CryoGrid (Westermann et al. 2022), or more specifically its sibling Julia implementation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently this package provides implementations of the following freeze curves:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Name Description Independent variable\nFreeWater Simple, piecewise-linear, enthalpy-based freeze/thaw scheme for pure \"free\" water. Enthalpy\nPainterKarra Coupled temperature-water retention model of Painter et al. (2014) Temperature\nDallAmico Coupled temperature-water retention model of Dall'Amico (2011). Temperature\nDallAmicoSalt Same as DallAmico but accounting for freezing point depressions due to salinity. Temperature\nMcKenzie Exponential-type empirical model of McKenzie et al. (2007) Temperature\nWestermann Nonlinear empirical model of Westermann et al. (2011) Temperature","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FreezeCurves.jl can be installed using the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add FreezeCurves","category":"page"},{"location":"","page":"Home","title":"Home","text":"or","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(\"FreezeCurves\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"All freezing characteristic curves are implemented as \"callable\" structs subtyping SFCCFunction. Callable structs are those with corresponding method definitions that allow them to be invoked like a function. As an example, we can initialize the freeze curve of McKenzie et al. (2007) with default parameter settings:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using FreezeCurves\n\nmkfc = McKenzie()","category":"page"},{"location":"","page":"Home","title":"Home","text":"The curve can then be evaluated at a particular temperature by calling it as a function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc(-1.0u\"°C\")\n# output\n1.860037988010418e-44","category":"page"},{"location":"","page":"Home","title":"Home","text":"The McKenzie struct contains default parameter settings necessary to evaluate the freeze curve. These can also be customized on construction. We could, for example, use a different value for the shape parameter, γ [K]","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc2 = McKenzie(γ=1.0u\"K\")\nmkfc2(-1.0u\"°C\")\n# output\n0.18393972058572117","category":"page"},{"location":"","page":"Home","title":"Home","text":"which results in a much longer curve with more than 18% unfrozen water content by volume at -1.0 °C!","category":"page"},{"location":"","page":"Home","title":"Home","text":"All freeze curves also accept their parameters (excluding some constants) as function arguments, with default values set to those found in the struct. For example, the full method signature for McKenzie is:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(f::McKenzie)(T; θtot, θsat, θres, Tₘ, γ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where θtot is the total (volumetric) water content (water + ice), θsat is the saturated water content (a.k.a natural porosity), θres is the residual water content, Tₘ is the melting temperature, and γ is the aforementioned shape parameter. We could have, therefore, adjusted the shape parameter by simply setting the keyword argument for our original mkfc:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc(-1.0u\"°C\", γ=1.0u\"K\")\n# output\n0.18393972058572117","category":"page"},{"location":"","page":"Home","title":"Home","text":"In cases where it is more convenient/appropriate to modify the default parameter values in the struct itself, it's worth noting that some parameters are nested one or two levels deep, e.g. setting the saturated water content (porosity) for a DallAmico freeze curve would look like:","category":"page"},{"location":"","page":"Home","title":"Home","text":"dafc = DallAmico(swrc=VanGenuchten(water=SoilWaterVolume(θsat=0.75)))","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is admittedly somewhat unwieldy, so FreezeCurves re-exports @set! from Setfield.jl to make this a bit easier:","category":"page"},{"location":"","page":"Home","title":"Home","text":"dafc = DallAmico()\n@set! dafc.swrc.vol.θsat = 0.75\ndafc.swrc.vol.θsat\n# output\n0.75","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note also that FreezeCurves.jl strictly uses unitful quantities (from Unitful.jl) to ensure physical coherence and avoid confusion between temperature units. The macro u\"°C\" applies units of degrees Celsius to the floating point number -0.1 (the degree symbol can be typed in the Julia REPL or editor using \\degree followed the TAB key). Speical method dispatches for ustrip are provided for freeze and water retention curve function types that will automatically strip the units from all fields nested within the function struct. When used without units, all temperatures are assumed to be in degrees Celsius (°C)!","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc_nounits = ustrip(mkfc)\nmkfc_nounits(-0.1)\n# output\n0.1839397205856375","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Brooks RH. Hydraulic properties of porous media. Colorado State University; 1965.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Dall'Amico M, Endrizzi S, Gruber S, Rigon RJ. A robust and energy-conserving model of freezing variably-saturated soil. The Cryosphere. 2011 Jun 1;5(2):469-84.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Koopmans RW, Miller RD. Soil freezing and soil water characteristic curves. Soil Science Society of America Journal. 1966 Nov;30(6):680-5.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Kurylyk BL, Watanabe K. The mathematical representation of freezing and thawing processes in variably-saturated, non-deformable soils. Advances in Water Resources. 2013 Oct 1;60:160-77.","category":"page"},{"location":"","page":"Home","title":"Home","text":"McKenzie JM, Voss CI, Siegel DI. Groundwater flow with energy transport and water–ice phase change: numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in water resources. 2007 Apr 1;30(4):966-83.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Westermann S, Boike J, Langer M, Schuler TV, Etzelmüller B. Modeling the impact of wintertime rain events on the thermal regime of permafrost. The Cryosphere. 2011 Oct 26;5(4):945-59.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Westermann S, Ingeman-Nielsen T, Scheer J, Aalstad K, Aga J, Chaudhary N, Etzelmüller B, Filhol S, Kääb A, Renette C, Schmidt LS. The CryoGrid community model (version 1.0)–a multi-physics toolbox for climate-driven simulations in the terrestrial cryosphere. Geoscientific Model Development Discussions. 2022 Jun 7:1-61.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Van Genuchten MT. A closed‐form equation for predicting the hydraulic conductivity of unsaturated soils. Soil science society of America journal. 1980 Sep;44(5):892-8.","category":"page"},{"location":"inference/#Estimating-freeze-curve-parameters","page":"Inference","title":"Estimating freeze curve parameters","text":"","category":"section"},{"location":"inference/","page":"Inference","title":"Inference","text":"Most soil freezing characteristic curves have one or more parameters which determine the shape of the curve. These parameters are typically determined empirically and can vary considerably between different types of soil.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"FreezeCurves.jl facilitates inference of these parameters via the Turing.jl probabilistic programming language. In order to keep FreezeCurves.jl as lightweight as possible, Turing is not included as a dependncy by default. It must be installed separately (Pkg.add(\"Turing\")) to expose the FreezeCurves.Inference module.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"We can demonstrate this with an idealized test case where simply corrupt the \"true\" liquid water content values with isotropic Gaussian noise. Here we fit the van Genuchten parameters for the Dall'Amico (2011) freezing characteristic curve using the Bayesian probabilistic freeze curve model pre-sepcified in the Inference module. This model assumes that the measurement errors are normally distributed, which matches our idealized test case.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"using Turing\nusing FreezeCurves\nusing FreezeCurves.Inference\n\nimport Random\n\nrng = Random.MersenneTwister(1234)\n\nfc = DallAmico(swrc=VanGenuchten(α=0.1u\"1/m\", n=1.8))\nTrange = vcat(-5.0u\"°C\":0.1u\"K\":-0.11u\"°C\", -0.1u\"°C\":0.001u\"K\":0.0u\"°C\")\nθtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))\nsfcc_model = SFCCModel(fc) # from FreezeCurves.Inference\nm = sfcc_model(Trange, θtrue) # condition on data\n# draw 1,000 samples using the No U-Turn Sampler; gradients are computed automatically by Turing using forward-mode automatic differentiation (ForwardDiff.jl).\nchain = sample(rng, m, NUTS(), 1000)\ndisplay(chain)","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Output:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Chains MCMC chain (1000×18×1 Array{Float64, 3}):\n\nIterations        = 501:1:1500\nNumber of chains  = 1\nSamples per chain = 1000\nWall duration     = 5.07 seconds\nCompute duration  = 5.07 seconds\nparameters        = logα, logn, Tₘ, θfrac, θres, σ\ninternals         = lp, n_steps, is_accept, acceptance_rate, log_density, hamiltonian_energy, hamiltonian_energy_error, max_hamiltonian_energy_error, tree_depth, numerical_error, step_size, nom_step_size\n\nSummary Statistics\n  parameters      mean       std   naive_se      mcse        ess      rhat   ess_per_sec \n      Symbol   Float64   Float64    Float64   Float64    Float64   Float64       Float64 \n\n        logα   -2.2728    0.0630     0.0020    0.0060    67.7835    1.0039       13.3590\n        logn   -0.2151    0.0513     0.0016    0.0039   145.5487    1.0033       28.6852\n          Tₘ   -0.0051    0.0036     0.0001    0.0004    45.3223    1.0085        8.9323\n       θfrac    0.4869    0.0049     0.0002    0.0004    80.6756    1.0039       15.8998\n        θres    0.0040    0.0032     0.0001    0.0002   326.9962    1.0007       64.4455\n           σ    0.0187    0.0012     0.0000    0.0001   243.4807    1.0007       47.9859\n\nQuantiles\n  parameters      2.5%     25.0%     50.0%     75.0%     97.5% \n      Symbol   Float64   Float64   Float64   Float64   Float64 \n\n        logα   -2.3786   -2.3181   -2.2779   -2.2364   -2.1379\n        logn   -0.3026   -0.2513   -0.2180   -0.1816   -0.1120\n          Tₘ   -0.0145   -0.0074   -0.0046   -0.0023   -0.0002\n       θfrac    0.4778    0.4836    0.4870    0.4901    0.4962\n        θres    0.0001    0.0014    0.0034    0.0060    0.0117\n           σ    0.0167    0.0179    0.0186    0.0195    0.0214","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"This implementation uses the log-transform of the van Genuchten parameters which, while not strictly necessary, has a few positive effects:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"It transforms the shape parameters, which have bounded support, into unconstrained space allowing for simple, unconstrained Gaussian priors.\nIt linearizes the otherwise non-linear relationship between α and n.\nIt smooths the partial derivatives of the objective function w.r.t the parameters which can improve the performance of gradient-based sampling and optimization methods (NUTS is a variant of Hamiltonian Monte Carlo which uses the gradient).","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"We can manually obtain the mean estimates of original parameters as α = mean(exp.(Array(group(chain, :logα)))) ≈ 0.103 and n = mean(1 .+ exp.(Array(group(chain, :logn)))) ≈ 1.808, which are quite close to the \"true\" values as expected.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Real measurement data can be used by simply replacing Trange and θtrue in this example with equal-length vectors of temperature and water content measurements. In cases where the isotropic Gaussian error assumption does not hold, the SFCCModel interface can be extended to use custom model implementations.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"If you aren't interested in the full posterior distribution, Turing also provides a very convenient interface to get a maximum a posteriori (MAP) estimate (i.e. a mode of the posterior) using Optim.jl:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"using Optim\n\noptimize(m, MAP(), LBFGS())","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Output:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"ModeResult with maximized lp of 394.29\n6-element Named Vector{Float64}\nA      │ \n───────┼────────────\n:logα  │    -2.28583\n:logn  │    -0.24065\n:Tₘ    │ -0.00340876\n:θfrac │    0.490879\n:θres  │ 2.05868e-21\n:σ     │    0.018277","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Alternatively, one can ignore the prior entirely and just get a maximum likelihood estimate:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"# here we use the common LBFGS optimizer; see Optim.jl docs for more options\noptimize(m, MLE(), LBFGS())","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Output:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"ModeResult with maximized lp of 394.29\n6-element Named Vector{Float64}\nA      │ \n───────┼────────────\n:logα  │    -2.28583\n:logn  │    -0.24065\n:Tₘ    │ -0.00340876\n:θfrac │    0.490879\n:θres  │ 8.58955e-23\n:σ     │    0.018277","category":"page"}]
}
