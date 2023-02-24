"""
    plot(T, fc::SFCC)

Simple plotting recipe for `SFCC` which takes care of stripping units.
"""
@recipe function plot(T, fc::SFCC)
    x := ustrip.(u"°C", T)
    xlabel := "Temperature (°C)"
    y := map(ustrip(fc) ∘ Base.Fix1(ustrip, u"°C"), T)
    ylabel := "Volumetric liquid water content"
    leg := :topleft
    ()
end
"""
    plot(Ts, fc::SFCC)

Simple plotting recipe for `SWRC` which takes care of stripping units.
"""
@recipe function plot(ψ, swrc::SWRC)
    x := ustrip.(u"m", ψ)
    xlabel := "Water matric potential (m)"
    y := map(ustrip(swrc) ∘ Base.Fix1(ustrip, u"m"), ψ)
    ylabel := "Volumetric liquid water content"
    leg := :topleft
    ()
end