# Free water freeze curve
"""
    FreeWater <: FreezeCurve

"Free water" freeze curve in terms of enthalpy (H), total water content (θtot), and
the latent heat of fusion of water (L).
"""
struct FreeWater <: FreezeCurve end

# for backwards compatibility
freewater(H, θtot, L) = freewater(H; θtot, L)
function freewater(H; θtot=0.5, θres=0.0, L=3.34e8)
    θtot = max(1e-8, θtot - θres)
    Lθ = L*θtot
    θw = IfElse.ifelse(
        H < zero(θtot),
        # Case 1: H < 0 -> frozen
        θres,
        # Case 2: H >= 0
        IfElse.ifelse(
            H >= Lθ,
            # Case 2a: H >= Lθ -> thawed
            θtot,
            # Case 2b: 0 <= H < Lθ -> phase change
            θres + H / Lθ
        )
    )
    return (;θw, Lθ)
end
(freeW::FreeWater)(H; θtot, θres=0.0, L=3.34e8) = freewater(H; θtot, θres, L).θw
