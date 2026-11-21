# Scalar test cases `name => (f, x₀)` for differentiability tests, where `f` maps a single input or parameter to
# the output of a freeze or water retention curve. All curves are stripped of units and temperatures are in °C.

# central finite difference used as the reference derivative
finitediff(f, x; h=1e-6) = (f(x + h) - f(x - h)) / (2h)

function autodiff_cases()
    vg = ustrip(VanGenuchten(α=2.0u"1/m", n=1.5))
    bc = ustrip(BrooksCorey(ψₛ=0.1u"m", λ=0.5))
    vol = (θsat=0.4, θres=0.05)
    ψ₀ = -1.0
    T₀ = -0.5
    swrc_cases = [
        "VanGenuchten: ψ" => (ψ -> vg(ψ; vol...), ψ₀),
        "VanGenuchten: α" => (α -> vg(ψ₀; vol..., α), 2.0),
        "VanGenuchten: n" => (n -> vg(ψ₀; vol..., n), 1.5),
        "VanGenuchten: θsat" => (θsat -> vg(ψ₀; θsat, θres=0.05), 0.4),
        "VanGenuchten: θres" => (θres -> vg(ψ₀; θsat=0.4, θres), 0.05),
        "VanGenuchten inverse: θ" => (θ -> inv(vg)(θ; vol...), 0.2),
        "BrooksCorey: ψ" => (ψ -> bc(ψ; vol...), ψ₀),
        "BrooksCorey: ψₛ" => (ψₛ -> bc(ψ₀; vol..., ψₛ), 0.1),
        "BrooksCorey: λ" => (λ -> bc(ψ₀; vol..., λ), 0.5),
        "BrooksCorey inverse: θ" => (θ -> inv(bc)(θ; vol...), 0.2),
    ]
    sfcc_cases = [
        "McKenzie: T" => (T -> ustrip(McKenzie())(T; vol..., γ=0.5), T₀),
        "McKenzie: γ" => (γ -> ustrip(McKenzie())(T₀; vol..., γ), 0.5),
        "Westermann: T" => (T -> ustrip(Westermann())(T; vol...), T₀),
        "Westermann: δ" => (δ -> ustrip(Westermann())(T₀; vol..., δ), 0.1),
        "Langer: T" => (T -> ustrip(Langer())(T; vol...), T₀),
        "Langer: a" => (a -> ustrip(Langer())(T₀; vol..., a), 19.0),
        "Hu2020: T" => (T -> ustrip(Hu2020())(T; vol...), T₀),
        "Hu2020: b" => (b -> ustrip(Hu2020())(T₀; vol..., b), 0.01),
        "PowerLaw: T" => (T -> ustrip(PowerLaw())(T; vol...), T₀),
        "PowerLaw: β" => (β -> ustrip(PowerLaw())(T₀; vol..., β), 0.5),
        "DallAmico: T" => (T -> ustrip(DallAmico())(T; vol...), T₀),
        "DallAmico: sat" => (sat -> ustrip(DallAmico())(0.5, sat; vol...), 0.8),
        "DallAmico: α" => (α -> ustrip(DallAmico())(T₀; vol..., α, n=1.5), 2.0),
        "DallAmico: n" => (n -> ustrip(DallAmico())(T₀; vol..., α=2.0, n), 1.5),
        "DallAmico: Tₘ" => (Tₘ -> ustrip(DallAmico())(T₀; vol..., Tₘ), -0.1),
        "PainterKarra: T" => (T -> ustrip(PainterKarra())(T; vol...), T₀),
        "PainterKarra: β" => (β -> ustrip(PainterKarra())(T₀; vol..., β, ω=0.5), 1.0),
        "PainterKarra: ω" => (ω -> ustrip(PainterKarra())(T₀; vol..., ω), 0.5),
        "DallAmicoSalt: T" => (T -> ustrip(DallAmicoSalt())(T; vol...), -2.0),
        "DallAmicoSalt: saltconc" => (saltconc -> ustrip(DallAmicoSalt())(-2.0; vol..., saltconc), 500.0),
    ]
    # differentiability of the analytical SFCC derivatives ∂θw∂T themselves
    freezecurves = (McKenzie(γ=0.5u"K"), Westermann(), Langer(), Hu2020(), PowerLaw(), DallAmico(), PainterKarra(β=1.0, ω=0.5))
    sfcc_derivative_cases = [
        "$(nameof(typeof(fc))) derivative: T" => (T -> FreezeCurves.derivative(ustrip(fc))(T; vol...), T₀) for fc in freezecurves
    ]
    push!(sfcc_derivative_cases, "DallAmicoSalt derivative: T" => (T -> FreezeCurves.derivative(ustrip(DallAmicoSalt()))(T; vol...), -2.0))
    return vcat(swrc_cases, sfcc_cases, sfcc_derivative_cases)
end
