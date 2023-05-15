var documenterSearchIndex = {"docs":
[{"location":"api/Solvers/#Solvers","page":"Solvers","title":"Solvers","text":"","category":"section"},{"location":"api/Solvers/","page":"Solvers","title":"Solvers","text":"The Solvers sub-module provides an interface for performing non-linear root-finding, e.g. solving implicit freeze curve formulations or recovering temperature given enthalpy, on the soil freeze characteristic curves implemented in this package. Note that this is not to be confused with parameter inference which is handled by a separate sub-module.","category":"page"},{"location":"api/Solvers/","page":"Solvers","title":"Solvers","text":"DocTestSetup = quote\n    using FreezeCurves\n    using FreezeCurves.Solvers\nend","category":"page"},{"location":"api/Solvers/","page":"Solvers","title":"Solvers","text":"Modules = [FreezeCurves.Solvers]\nPrivate = false\nOrder = [:type, :function, :macro]","category":"page"},{"location":"api/Solvers/#FreezeCurves.Solvers.LUT","page":"Solvers","title":"FreezeCurves.Solvers.LUT","text":"LUT{N,T,coordnames,TCoords,TData}\n\nGeneric implementation of an N-dimensional lookup table with support for reverse lookups w.r.t one variable.\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.SFCCInverseEnthalpyObjective","page":"Solvers","title":"FreezeCurves.Solvers.SFCCInverseEnthalpyObjective","text":"SFCCInverseEnthalpyObjective{TF,Tkwargs<:NamedTuple,Thc,TL,TH,Tsat} <: AbstractSFCCObjective\n\nOptimization objective for finding a temperature, T, wich resolves the conservation equation:\n\n0 = (H - Lθ)C - T\n\ngiven a fixed enthalpy value H and freeze curve arguments f_args.\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.SFCCNewtonSolver","page":"Solvers","title":"FreezeCurves.Solvers.SFCCNewtonSolver","text":"SFCCNewtonSolver <: SFCCSolver\n\nFast, specialized implementation of Newton's method with backtracking line search for resolving the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid jumping over the solution. This prevents convergence issues that arise due to discontinuities and strong non-linearity in most common soil freeze curves.\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.SFCCPreSolver","page":"Solvers","title":"FreezeCurves.Solvers.SFCCPreSolver","text":"SFCCPreSolver{TCache} <: SFCCSolver\n\nA fast SFCC \"solver\" which pre-builds an interpolant for the freeze curve in terms of enthalpy, θ(H).\n\n\n\n\n\n","category":"type"},{"location":"api/Solvers/#FreezeCurves.Solvers.build_lut-Tuple{Any, Vararg{Pair{Symbol, var\"#s48\"} where var\"#s48\"<:Union{Number, AbstractVector{T} where T}, N} where N}","page":"Solvers","title":"FreezeCurves.Solvers.build_lut","text":"build_lut(f, coords::Pair{Symbol,<:Union{Number,AbstractVector}}...; f_kwargs...)\n\nBuilds an n-dimensional lookup table LUT by invoking f on the given coords. f should be a function which accepts keyword arguments matching coords as well as any additional constant keyword arguments f_kwargs.\n\nExample:\n\nf(; x=1.0, y=1.0) = x + y\nlut = build_lut(f, :x => -10.0:10.0, :y=-10.0:10.0)\nlut(; x=2.2, y=1.5)\n# output\n3.7\n\n\n\n\n\n","category":"method"},{"location":"api/Solvers/#FreezeCurves.Solvers.sfccsolve-Union{Tuple{return_all}, Tuple{FreezeCurves.Solvers.AbstractSFCCObjective, SFCCSolver, Any}, Tuple{FreezeCurves.Solvers.AbstractSFCCObjective, SFCCSolver, Any, Val{return_all}}} where return_all","page":"Solvers","title":"FreezeCurves.Solvers.sfccsolve","text":"sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀, ::Val{return_all}=Val{true}()) where {return_all}\n\nSolve the given objective obj using solver and initial guess x₀. If return_all=true, then sfccsolve should return a named tuple with at least the temperature solution T, the liquid water content θw, the heat capacity C, and the liquid water content deriviative ∂θw∂T defined. Solver-specific additional values may also be included.\n\n\n\n\n\n","category":"method"},{"location":"api/Solvers/#FreezeCurves.Solvers.sfccsolve-Union{Tuple{return_all}, Tuple{FreezeCurves.Solvers.SFCCInverseEnthalpyObjective, FreezeCurves.Solvers.SFCCNewtonSolver, Number}, Tuple{FreezeCurves.Solvers.SFCCInverseEnthalpyObjective, FreezeCurves.Solvers.SFCCNewtonSolver, Number, Val{return_all}}} where return_all","page":"Solvers","title":"FreezeCurves.Solvers.sfccsolve","text":"sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNewtonSolver, T₀::Number, ::Val{return_all}=Val{true}(); error_when_not_converged=true)\n\nSolves obj using the specialized Newton solver and returns the result. If return_all is true, a named tuple (;T, Tres, θw, C, itercount) is returned; otherwise (by default), only the temperature solution is returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#Index-of-public-API","page":"Index of public API","title":"Index of public API","text":"","category":"section"},{"location":"api/","page":"Index of public API","title":"Index of public API","text":"","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"api/FreezeCurves/#FreezeCurves","page":"FreezeCurves","title":"FreezeCurves","text":"","category":"section"},{"location":"api/FreezeCurves/","page":"FreezeCurves","title":"FreezeCurves","text":"Modules = [FreezeCurves]\nPrivate = false\nOrder = [:type, :function, :macro]","category":"page"},{"location":"api/FreezeCurves/#FreezeCurves.BrooksCorey","page":"FreezeCurves","title":"FreezeCurves.BrooksCorey","text":"BrooksCorey{Twp,Tψₛ,Tλ} <: SWRC\n\nvan Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.     Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.DallAmico","page":"FreezeCurves","title":"FreezeCurves.DallAmico","text":"DallAmico{TFT,Tg,Tswrc<:SWRC} <: SFCC\n\nDall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.DallAmicoSalt","page":"FreezeCurves","title":"FreezeCurves.DallAmicoSalt","text":"DallAmicoSalt{TFT,Tsc,TR,Tg,Tswrc<:SWRC} <: SFCC\n\nFreeze curve from Dall'Amico (2011) with saline freezing point depression.\n\nAngelopoulos M, Westermann S, Overduin P, Faguet A, Olenchenko V, Grosse G, Grigoriev MN. Heat and salt flow in subsea permafrost     modeled with CryoGRID2. Journal of Geophysical Research: Earth Surface. 2019 Apr;124(4):920-37.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.FreeWater","page":"FreezeCurves","title":"FreezeCurves.FreeWater","text":"FreeWater <: FreezeCurve\n\n\"Free water\" freeze curve in terms of enthalpy (H), total water content (θtot), and the latent heat of fusion of water (L).\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.FreezeCurve","page":"FreezeCurves","title":"FreezeCurves.FreezeCurve","text":"FreezeCurve\n\nBase type for freezing characteristic curves which relate temperature or enthalpy to unfrozen water content.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.Hu2020","page":"FreezeCurves","title":"FreezeCurves.Hu2020","text":"Hu2020{TFT,Tvol,Tb}  <: SFCC\n\nSoil freezing characteristic curve formulation of Hu et al. 2020.\n\nHu G, Zhao L, Zhu X, Wu X, Wu T, Li R, Xie C, Hao J. Review of algorithms and parameterizations to determine unfrozen water content in frozen soil. Geoderma. 2020 Jun 1;368:114277.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.McKenzie","page":"FreezeCurves","title":"FreezeCurves.McKenzie","text":"McKenzie{TFT,Tvol,Tγ} <: SFCC\n\nMcKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:     numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,     30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.PainterKarra","page":"FreezeCurves","title":"FreezeCurves.PainterKarra","text":"PainterKarra{TFT,Tβ,Tω,Tg,Tswrc<:SWRC} <: SFCC\n\nPainter SL, Karra S. Constitutive Model for Unfrozen Water Content in Subfreezing Unsaturated Soils. Vadose zone j. 2014 Apr;13(4):1-8.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.PowerLaw","page":"FreezeCurves","title":"FreezeCurves.PowerLaw","text":"PowerLaw{Tvol,Tα,Tβ} <: SFCC\n\nCommonly used power law relation of Lovell (1957).\n\nLovell, C.: Temperature effects on phase composition and strength of partially frozen soil, Highway Research Board Bulletin, 168, 74–95, 1957.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SFCC","page":"FreezeCurves","title":"FreezeCurves.SFCC","text":"SFCC <: FreezeCurve\n\nBase type for a soil freeze characteristic curve (SFCC) function. Subtypes should be callable structs that implement the freeze curve and contain any necessary additional constants or configuration options. User-specified parameters can either be supplied in the struct or declared as model parameters via the variables method.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SFCCSolver","page":"FreezeCurves","title":"FreezeCurves.SFCCSolver","text":"SFCCSolver\n\nBase type representing non-linear solvers for implicit SFCC functions.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SWRC","page":"FreezeCurves","title":"FreezeCurves.SWRC","text":"SWRC\n\nBase type for soil water retention curves (SWRC) which relate soil water matric potential to water content.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SoilFreezeThawProperties","page":"FreezeCurves","title":"FreezeCurves.SoilFreezeThawProperties","text":"SoilFreezeThawProperties{TTₘ,TLsl}\n\nStruct containing constants/parameters common to some or all SFCCs.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SoilFreezeThawProperties-Tuple{SFCC}","page":"FreezeCurves","title":"FreezeCurves.SoilFreezeThawProperties","text":"SoilFreezeThawProperties(f::SFCC)\n\nRetrieves the default SoilFreezeThawProperties from f; should be defined for all freeze curves.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.SoilWaterVolume","page":"FreezeCurves","title":"FreezeCurves.SoilWaterVolume","text":"SoilWaterVolume{Tρw,Tθres,Tθsat}\n\nStruct containing basic physical constants and propeties related to soil pore water.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.SoilWaterVolume-Tuple{SFCC}","page":"FreezeCurves","title":"FreezeCurves.SoilWaterVolume","text":"SoilWaterVolume(f::SFCC)\n\nRetrieves the nested SoilWaterVolume from the SoilFreezeThawProperties of the freeze curve f.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.SoilWaterVolume-Tuple{SWRC}","page":"FreezeCurves","title":"FreezeCurves.SoilWaterVolume","text":"SoilWaterVolume(f::SWRC)\n\nGet the SoilWaterVolume defined for the given SWRC f. Must be implemented for all SWRC types. Default implementation retrieves the field water.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.VanGenuchten","page":"FreezeCurves","title":"FreezeCurves.VanGenuchten","text":"VanGenuchten{Tvol,Tα,Tn} <: SWRC\n\nvan Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.     Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.\n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.Westermann","page":"FreezeCurves","title":"FreezeCurves.Westermann","text":"Westermann{TFT,Tvol,Tδ} <: SFCC\n\nWestermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of     wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,     https://doi.org/10.5194/tc-5-945-2011, 2011. \n\n\n\n\n\n","category":"type"},{"location":"api/FreezeCurves/#FreezeCurves.inflectionpoint-Tuple{SFCC}","page":"FreezeCurves","title":"FreezeCurves.inflectionpoint","text":"inflectionpoint(f::SFCC)\n\nReturns the analytical inflection point (i.e. where ∂²θ/∂T^2 = 0), if available.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#FreezeCurves.inflectionpoint-Tuple{SWRC}","page":"FreezeCurves","title":"FreezeCurves.inflectionpoint","text":"inflectionpoint(f::SWRC)\n\nReturns the analytical solution for the inflection point (i.e. where ∂²θ/∂ψ² = 0) of the SWRC, if available.\n\n\n\n\n\n","category":"method"},{"location":"api/FreezeCurves/#Unitful.ustrip-Tuple{Union{SFCC, SWRC}}","page":"FreezeCurves","title":"Unitful.ustrip","text":"ustrip(x)\n\nReconstructs the type or function x with all numerical quantities stripped of units.\n\n\n\n\n\n","category":"method"},{"location":"#FreezeCurves.jl","page":"Home","title":"FreezeCurves.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FreezeCurves is a lightweight Julia package to facilitate the study and analysis of soil freezing characteristics (Koopmans and Miller, 1966). The relationship between temperature and unfrozen water content in porous media (such as soils) is often highly nonlinear and plays a significant role in the analysis of freeze-thaw dynamics in science and engineering. One common application is in the geophysical modeling of permafrost, where having faithful representations of freeze-thaw processes is often paramount to accurately resolving long-term changes in the subsurface thermal regime.","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"index.md\",\"inference.md\",\"api/index.md\"]","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The soil freezing characteristic curve (SFCC) is typically a monotonic function which maps (usually subzero) temperatures [°C] to volumetric unfrozen/liquid water contents [m³/m³]. The SFCC is closely related to the soil-water retention curve (SWRC) which governs the similar non-linear realtionship between soil-water matric potential [m] saturation [m³/m³]. SWRCs are widely used in hydrological modeling, with the two most common formulations being the van Genuchten (1980) and Brooks-Corey (Brooks, 1965) models. SFCCs are less widely known, but numerous models have been proposed over the years (Kurylyk and Watanbe, 2013).","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package is intended to act as a living repository of SFCC/SWRC implementations which can then be fitted to data or consumed downstream by thermal or hydrological models such as CryoGrid (Westermann et al. 2022), or more specifically its sibling Julia implementation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently this package provides implementations of the following freeze curves:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Name Description Independent variable\nFreeWater Simple, piecewise-linear, enthalpy-based freeze/thaw scheme for pure \"free\" water. Enthalpy\nPainterKarra Coupled temperature-water retention model of Painter et al. (2014) Temperature\nDallAmico Coupled temperature-water retention model of Dall'Amico (2011). Temperature\nDallAmicoSalt Same as DallAmico but accounting for freezing point depressions due to salinity. Temperature\nMcKenzie Exponential-type empirical model of McKenzie et al. (2007) Temperature\nWestermann Nonlinear empirical model of Westermann et al. (2011) Temperature","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FreezeCurves.jl can be installed using the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add FreezeCurves","category":"page"},{"location":"","page":"Home","title":"Home","text":"or","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(\"FreezeCurves\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"All freezing characteristic curves are implemented as \"callable\" structs subtyping SFCC. Callable structs are those with corresponding method definitions that allow them to be invoked like a function. As an example, we can initialize the freeze curve of McKenzie et al. (2007) with default parameter settings:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using FreezeCurves\n\nmkfc = McKenzie()","category":"page"},{"location":"","page":"Home","title":"Home","text":"The curve can then be evaluated at a particular temperature by calling it as a function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc(-1.0u\"°C\")\n# output\n1.860037988010418e-44","category":"page"},{"location":"","page":"Home","title":"Home","text":"The McKenzie struct contains default parameter settings necessary to evaluate the freeze curve. These can also be customized on construction. We could, for example, use a different value for the shape parameter, γ [K]","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc2 = McKenzie(γ=1.0u\"K\")\nmkfc2(-1.0u\"°C\")\n# output\n0.18393972058572117","category":"page"},{"location":"","page":"Home","title":"Home","text":"which results in a much longer curve with more than 18% unfrozen water content by volume at -1.0 °C!","category":"page"},{"location":"","page":"Home","title":"Home","text":"All freeze curves also accept their parameters (excluding some constants) as function arguments, with default values set to those found in the struct. For example, the full method signature for McKenzie is:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(f::McKenzie)(T; θtot, θsat, θres, Tₘ, γ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where θtot is the total (volumetric) water content (water + ice), θsat is the saturated water content (a.k.a natural porosity), θres is the residual water content, Tₘ is the melting temperature, and γ is the aforementioned shape parameter. We could have, therefore, adjusted the shape parameter by simply setting the keyword argument for our original mkfc:","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc(-1.0u\"°C\", γ=1.0u\"K\")\n# output\n0.18393972058572117","category":"page"},{"location":"","page":"Home","title":"Home","text":"In cases where it is more convenient/appropriate to modify the default parameter values in the struct itself, it's worth noting that some parameters are nested one or two levels deep, e.g. setting the saturated water content (porosity) for a DallAmico freeze curve would look like:","category":"page"},{"location":"","page":"Home","title":"Home","text":"dafc = DallAmico(swrc=VanGenuchten(water=SoilWaterVolume(θsat=0.75)))","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is admittedly somewhat unwieldy, so FreezeCurves re-exports @set! from Setfield.jl to make this a bit easier:","category":"page"},{"location":"","page":"Home","title":"Home","text":"dafc = DallAmico()\n@set! dafc.swrc.vol.θsat = 0.75\ndafc.swrc.vol.θsat\n# output\n0.75","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note also that FreezeCurves.jl strictly uses unitful quantities (from Unitful.jl) to ensure physical coherence and avoid confusion between temperature units. The macro u\"°C\" applies units of degrees Celsius to the floating point number -0.1 (the degree symbol can be typed in the Julia REPL or editor using \\degree followed the TAB key). Speical method dispatches for ustrip are provided for freeze and water retention curve function types that will automatically strip the units from all fields nested within the function struct. When used without units, all temperatures are assumed to be in degrees Celsius (°C)!","category":"page"},{"location":"","page":"Home","title":"Home","text":"mkfc_nounits = ustrip(mkfc)\nmkfc_nounits(-0.1)\n# output\n0.1839397205856375","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Brooks RH. Hydraulic properties of porous media. Colorado State University; 1965.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Dall'Amico M, Endrizzi S, Gruber S, Rigon RJ. A robust and energy-conserving model of freezing variably-saturated soil. The Cryosphere. 2011 Jun 1;5(2):469-84.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Koopmans RW, Miller RD. Soil freezing and soil water characteristic curves. Soil Science Society of America Journal. 1966 Nov;30(6):680-5.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Kurylyk BL, Watanabe K. The mathematical representation of freezing and thawing processes in variably-saturated, non-deformable soils. Advances in Water Resources. 2013 Oct 1;60:160-77.","category":"page"},{"location":"","page":"Home","title":"Home","text":"McKenzie JM, Voss CI, Siegel DI. Groundwater flow with energy transport and water–ice phase change: numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in water resources. 2007 Apr 1;30(4):966-83.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Westermann S, Boike J, Langer M, Schuler TV, Etzelmüller B. Modeling the impact of wintertime rain events on the thermal regime of permafrost. The Cryosphere. 2011 Oct 26;5(4):945-59.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Westermann S, Ingeman-Nielsen T, Scheer J, Aalstad K, Aga J, Chaudhary N, Etzelmüller B, Filhol S, Kääb A, Renette C, Schmidt LS. The CryoGrid community model (version 1.0)–a multi-physics toolbox for climate-driven simulations in the terrestrial cryosphere. Geoscientific Model Development Discussions. 2022 Jun 7:1-61.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Van Genuchten MT. A closed‐form equation for predicting the hydraulic conductivity of unsaturated soils. Soil science society of America journal. 1980 Sep;44(5):892-8.","category":"page"},{"location":"inference/#Estimating-freeze-curve-parameters","page":"Inference","title":"Estimating freeze curve parameters","text":"","category":"section"},{"location":"inference/","page":"Inference","title":"Inference","text":"Most soil freezing characteristic curves have one or more parameters which determine the shape of the curve. These parameters are typically determined empirically and can vary considerably between different types of soil.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"FreezeCurves.jl facilitates inference of these parameters via the Turing.jl probabilistic programming language. In order to keep FreezeCurves.jl as lightweight as possible, Turing is not included as a dependncy by default. It must be installed separately (Pkg.add(\"Turing\")) to expose the FreezeCurves.Inference module.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"We can demonstrate this with an idealized test case where simply corrupt the \"true\" liquid water content values with isotropic Gaussian noise. Here we fit the van Genuchten parameters for the Dall'Amico (2011) freezing characteristic curve using the Bayesian probabilistic freeze curve model pre-sepcified in the Inference module. This model assumes that the measurement errors are normally distributed, which matches our idealized test case.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"using Turing\nusing FreezeCurves\nusing FreezeCurves.Inference\n\nimport Random\n\nrng = Random.MersenneTwister(1234)\n\nfc = DallAmico(swrc=VanGenuchten(α=0.1u\"1/m\", n=1.8))\nTrange = vcat(-5.0u\"°C\":0.1u\"K\":-0.11u\"°C\", -0.1u\"°C\":0.001u\"K\":0.0u\"°C\")\nθtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))\nsfcc_model = SFCCModel(fc) # from FreezeCurves.Inference\nm = sfcc_model(Trange, θtrue) # condition on data\n# draw 1,000 samples using the No U-Turn Sampler w/ 500 adaptation steps and 85% target acceptance rate; gradients are computed automatically by Turing using forward-mode automatic differentiation (ForwardDiff.jl).\nchain = sample(rng, m, NUTS(500,0.85), 1_000)\ndisplay(chain)","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Output:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Chains MCMC chain (1000×19×1 Array{Float64, 3}):\n\nIterations        = 501:1:1500\nNumber of chains  = 1\nSamples per chain = 1000\nWall duration     = 10.09 seconds\nCompute duration  = 10.09 seconds\nparameters        = logα, logn, Tₘ, sat, por, res, σ\ninternals         = lp, n_steps, is_accept, acceptance_rate, log_density, hamiltonian_energy, hamiltonian_energy_error, max_hamiltonian_energy_error, tree_depth, numerical_error, step_size, nom_step_size\n\nSummary Statistics\n  parameters      mean       std   naive_se      mcse        ess      rhat   ess_per_sec \n      Symbol   Float64   Float64    Float64   Float64    Float64   Float64       Float64 \n\n        logα   -2.2574    0.0575     0.0018    0.0033   381.2072    1.0010       37.7732\n        logn   -0.2221    0.0494     0.0016    0.0026   456.1224    0.9990       45.1964\n          Tₘ   -0.0034    0.0027     0.0001    0.0001   371.9088    0.9994       36.8518\n         sat    0.9889    0.0101     0.0003    0.0004   437.2432    1.0006       43.3257\n         por    0.4949    0.0052     0.0002    0.0002   368.9660    1.0042       36.5602\n         res    0.0074    0.0063     0.0002    0.0002   631.1171    1.0044       62.5364\n           σ    0.0193    0.0012     0.0000    0.0000   787.7707    0.9992       78.0589\n\nQuantiles\n  parameters      2.5%     25.0%     50.0%     75.0%     97.5% \n      Symbol   Float64   Float64   Float64   Float64   Float64 \n\n        logα   -2.3575   -2.2960   -2.2622   -2.2221   -2.1284\n        logn   -0.3147   -0.2585   -0.2219   -0.1887   -0.1282\n          Tₘ   -0.0100   -0.0050   -0.0028   -0.0012   -0.0001\n         sat    0.9621    0.9848    0.9919    0.9962    0.9997\n         por    0.4851    0.4915    0.4947    0.4980    0.5053\n         res    0.0002    0.0025    0.0058    0.0107    0.0238\n           σ    0.0171    0.0184    0.0192    0.0200    0.0218","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"This implementation uses the log-transform of the van Genuchten parameters which, while not strictly necessary, has a few positive effects:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"It transforms the shape parameters, which have bounded support, into unconstrained space allowing for simple, unconstrained Gaussian priors.\nIt linearizes the otherwise non-linear relationship between α and n.\nIt smooths the partial derivatives of the objective function w.r.t the parameters which can improve the performance of gradient-based sampling and optimization methods (NUTS is a variant of Hamiltonian Monte Carlo which uses the gradient).","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"We can manually obtain the mean estimates of original parameters as α = mean(exp.(Array(group(chain, :logα)))) ≈ 0.105 and n = mean(1 .+ exp.(Array(group(chain, :logn)))) ≈ 1.802, which are quite close to the \"true\" values as expected.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Real measurement data can be used by simply replacing Trange and θtrue in this example with equal-length vectors of temperature and water content measurements. In cases where the isotropic Gaussian error assumption does not hold, the SFCCModel interface can be extended to use custom model implementations.","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"If you aren't interested in the full posterior distribution, Turing also provides a very convenient interface to get a maximum a posteriori (MAP) estimate (i.e. a mode of the posterior) using Optim.jl:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"using Optim\n\noptimize(m, MAP(), LBFGS())","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Output:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"ModeResult with maximized lp of 394.16\n7-element Named Vector{Float64}\nA     │ \n──────┼─────────────\n:logα │      -2.3055\n:logn │    -0.228226\n:Tₘ   │ -1.45969e-52\n:sat  │          1.0\n:por  │     0.496235\n:res  │  4.83361e-31\n:σ    │     0.018805","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Alternatively, one can ignore the prior entirely and just get a maximum likelihood estimate:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"# here we use the common LBFGS optimizer; see Optim.jl docs for more options\noptimize(m, MLE(), LBFGS())","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Output:","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"ModeResult with maximized lp of 394.16\n7-element Named Vector{Float64}\nA     │ \n──────┼─────────────\n:logα │      -2.3055\n:logn │    -0.228226\n:Tₘ   │ -6.85422e-72\n:sat  │          1.0\n:por  │     0.496235\n:res  │  4.64697e-43\n:σ    │     0.018805","category":"page"},{"location":"inference/","page":"Inference","title":"Inference","text":"Note how the MLE vs MAP results are almost identical in this case since we have a (relatively) large generated dataset and the default priors are only very weakly informative.","category":"page"}]
}
