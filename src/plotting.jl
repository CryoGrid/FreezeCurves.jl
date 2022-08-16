"""
    plot(T, fc::SFCCFunction)

Simple plotting recipe for `SFCCFunction` which takes care of stripping units.
"""
@recipe function plot(T, fc::SFCCFunction)
    x := ustrip.(u"°C", T)
    xlabel := "Temperature (°C)"
    y := map(ustrip(stripparams(fc)) ∘ Base.Fix1(ustrip, u"°C"), T)
    ylabel := "Volumetric liquid water content"
    leg := :topleft
    ()
end
"""
    plot(Ts, fc::SFCCFunction)

Simple plotting recipe for `SWRCFunction` which takes care of stripping units.
"""
@recipe function plot(ψ, swrc::SWRCFunction)
    x := ustrip.(u"m", ψ)
    xlabel := "Water matric potential (m)"
    y := map(ustrip(stripparams(swrc)) ∘ Base.Fix1(ustrip, u"m"), ψ)
    ylabel := "Volumetric liquid water content"
    leg := :topleft
    ()
end