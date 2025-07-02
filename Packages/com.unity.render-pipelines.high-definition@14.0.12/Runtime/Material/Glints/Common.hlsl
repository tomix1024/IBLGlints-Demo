#ifndef GLINT_COMMON_HLSL
#define GLINT_COMMON_HLSL

float4 SanitizePositiveFinite(float4 x)
{
    return float4(
        SanitizePositiveFinite(x.x),
        SanitizePositiveFinite(x.y),
        SanitizePositiveFinite(x.z),
        SanitizePositiveFinite(x.w)
    );
}

float2x4 SanitizePositiveFinite(float2x4 x)
{
    return float2x4(
        SanitizePositiveFinite(x[0]),
        SanitizePositiveFinite(x[1])
    );
}

float4x4 SanitizePositiveFinite(float4x4 x)
{
    return float4x4(
        SanitizePositiveFinite(x[0]),
        SanitizePositiveFinite(x[1]),
        SanitizePositiveFinite(x[2]),
        SanitizePositiveFinite(x[3])
    );
}

#endif // GLINT_COMMON_HLSL
