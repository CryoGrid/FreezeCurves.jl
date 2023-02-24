# Free water freeze curve
"""
    FreeWater <: FreezeCurve

"Free water" freeze curve in terms of enthalpy (H), total water content (θtot), and
the latent heat of fusion of water (L).
"""
struct FreeWater <: FreezeCurve end
function freewater(H, θtot, L)
    θtot = max(1e-8, θtot)
    Lθ = L*θtot
    I_t = H > Lθ
    I_f = H <= 0.0
    I_c = (H > 0.0) && (H <= Lθ)
    # compute liquid water content -> heat capacity -> temperature
    θw = (I_c*(H/Lθ) + I_t)θtot
    return (;θw, I_t, I_f, I_c, Lθ)
end
(freeW::FreeWater)(H, θtot, L) = freewater(H, θtot, L).θw
