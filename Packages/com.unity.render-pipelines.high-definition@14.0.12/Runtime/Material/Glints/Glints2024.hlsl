
// Include Deliot & Belcour implementation
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Glints2023.hlsl"

float _Glintiness;

// Use uniform variabels for debugging...
float _ZeroIfPgt1;
//float _LogSinSunAngle;
float _SunSolidAngle;
int _GlintNDFIntegrationMode;

int _GlintVisualizeMultinomialWithLi;
int _GlintVisualizeMultinomialIndex;
int _GlintVisualizeMultinomialEx;

float ComputeTotalNDF(float roughness)
{
#if 1
    return 1.0 + pow(roughness, 0.5*2.821827587526826);
#else
    // https://www.wolframalpha.com/input?i=integrate+2*+a+%2F+%28%28x*a+-+x%29*x+%2B+1%29**2+dx+from+0+to+1
    // https://proofwiki.org/wiki/Arctangent_of_Imaginary_Number
    float a2 = Sq(roughness);
    float sqrt_term = sqrt(1-a2);
    // return 1 + 0.5*a2*log((1+sqrt_term)/(1-sqrt_term))/sqrt_term;
    return 1 + a2 * log((1 + sqrt_term) / roughness) / sqrt_term;
#endif
}

//#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/GlintsMultinomial.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/GlintsMultinomial8.hlsl"



float SampleGlints2024NDF(float3 halfwayTS, float LdotH, float roughness, float2 uv, float2 duvdx, float2 duvdy)
{
    float Dtarget = D_GGX(halfwayTS.z, roughness);

    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);
    // Sun:
    // omega = 6.8e-5 sr.
    // sin(gamma) = 4.6e-3 == gamma
    //float sinGammaSq = exp(2*_LogSinSunAngle);
    //float cosGamma = sqrt(1 - sinGammaSq);
    //float Al = 2*PI * (1.0 - cosGamma);
    float Al = _SunSolidAngle;
    float Ah = Al / abs(4*LdotH);
    float integratedNDF = Dtarget * Ah;

    float successProb = integratedNDF / totalNDF;
    if (successProb > 1)
        successProb *= 1-_ZeroIfPgt1;

    // Explicitly simulate smooth surface if _LogMicrofacetDensity < 0
    if (_LogMicrofacetDensity <= 0)
        return 1;

    // NOTE: This includes D(NdotH) factor from smooth model...
    float D = SampleGlints2023NDF_Internal(halfwayTS, successProb, uv, duvdx, duvdy);
    //return Dtarget * (D / successProb);
    return D * totalNDF / Ah;
}


float SampleGlints2024NDF_Area(float3 halfwayTS, float LdotH, float roughness, float integratedNDF, float4x3 lightVerts, float2 uv, float2 duvdx, float2 duvdy)
{
    // integratedNDF = \int_light D(H)/(4*LdotH) dL

    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);

    if (integratedNDF > totalNDF)
        integratedNDF *= 1-_ZeroIfPgt1;
    float p = saturate(integratedNDF / totalNDF); // = R*Dtarget/Dmax

    // Skip computation below if probability is zero!
    if (p == 0)
        return 0;

    // Explicitly simulate smooth surface if _LogMicrofacetDensity < 0
    if (_LogMicrofacetDensity <= 0)
        return 1;

    // Division by microfacet count handled internally.
    float D = SampleGlints2023NDF_Internal(halfwayTS, p, uv, duvdx, duvdy) / p;
    return D;
}


float SampleGlints2024NDF_IBL(float4 levels, float4 probs, float roughness, float visibleNDF, float2 uv, float2 duvdx, float2 duvdy)
{
    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);
    float baseProb = saturate(visibleNDF / totalNDF);

    // TODO total NDF vs visible NDF vs integrated NDF!!!

    // NOTE: assume levels.x == 0!
    float4 successProbs = baseProb * probs;
    successProbs.x += (1-baseProb);

    float4 Ds = SampleGlintsMultinomialNDF_Internal(successProbs, uv, duvdx, duvdy);
#ifdef _GLINTS_VISUALIZE_MULTINOMIAL
    float result = Ds[clamp(_GlintVisualizeMultinomialIndex, 0, 3)];
    if (_GlintVisualizeMultinomialEx != 0)
    {
        result = successProbs[clamp(_GlintVisualizeMultinomialIndex, 0, 3)];
    }
    if (_GlintVisualizeMultinomialWithLi != 0)
    {
        result *= levels[clamp(_GlintVisualizeMultinomialIndex, 0, 3)];
    }
#else // _GLINTS_VISUALIZE_MULTINOMIAL
    float result = dot(levels.xyzw, Ds.xyzw) / (baseProb * dot(levels.xyzw, probs.xyzw));
#endif // _GLINTS_VISUALIZE_MULTINOMIAL
    return result;
}

float SampleGlints2024NDF_IBL(float2x4 levels, float2x4 probs, float roughness, float visibleNDF, float2 uv, float2 duvdx, float2 duvdy)
{
    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);
    float baseProb = saturate(visibleNDF / totalNDF);

    // TODO total NDF vs visible NDF vs integrated NDF!!!

    // NOTE: assume levels.x == 0!
    float2x4 successProbs = baseProb * probs;
    successProbs[0][0] += (1-baseProb);
    levels[0][0] = 0;

    float2x4 Ds = SampleGlintsMultinomial8NDF_Internal(successProbs, uv, duvdx, duvdy);

    //float Dsum = (dot(Ds[0], 1) + dot(Ds[1], 1));
    //float psum = (dot(successProbs[0], 1) + dot(successProbs[1], 1));
    //return Dsum / psum;

    float result = (dot(levels[0], Ds[0])+dot(levels[1], Ds[1])) / (dot(levels[0], successProbs[0]) + dot(levels[1], successProbs[1]));
    return result;
}

float SampleGlints2024NDF_IBL_Wang2020(float roughness, float visibleNDF, float2 uv, float2 duvdx, float2 duvdy)
{
    // Carefully looking at the equations in Wang2020, we came to the following conclusions:
    // - the spatial traversal can be replaced with Delito's grid for speedup.
    // - the prefiltered Lo is supposed to work identically to our implementation.
    // - individual glints are simply sampled with p=|\Omega_o|. We simply reuse the _SunSolidAngle for directiona lights at this point...

    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);
    float successProb = saturate(_SunSolidAngle * visibleNDF / totalNDF);

    // Division by microfacet count handled internally.
    float3 unusedH = 0.0.rrr; // There is no halfway vector to interpolate over. In the compiler we trust!
    float result = SampleGlints2023NDF_Internal(unusedH, successProb, uv, duvdx, duvdy) / successProb;
    return result;
}
