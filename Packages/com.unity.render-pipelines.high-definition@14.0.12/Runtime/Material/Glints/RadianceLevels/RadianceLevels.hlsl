#ifndef PREINTEGRATED_LIGHT_DIR_HLSL
#define PREINTEGRATED_LIGHT_DIR_HLSL

#define GLINT_LEVELS_SLICE_COUNT 3

TEXTURECUBE_ARRAY(_SkyTextureGlintLevelWeights);

cbuffer _GlintLevelsData {
    float4 _GlintLevels8_0;
    float4 _GlintLevels8_1;
    float4 _GlintLevels4;
    float4 _GlintLevelsDummy;
};


void GetPreIntegratedLightLevelsFromSky(float3 iblR, float iblPerceptualRoughness, out float4 L, out float4 p, int sliceIndex = 0)
{
    float iblMipLevel = PerceptualRoughnessToMipmapLevel(iblPerceptualRoughness);
    p = SAMPLE_TEXTURECUBE_ARRAY_LOD(_SkyTextureGlintLevelWeights, s_trilinear_clamp_sampler, iblR, sliceIndex*GLINT_LEVELS_SLICE_COUNT + 2, iblMipLevel);
    // NOTE: Normalization is required for low-precision texture formats only.
    float psum = dot(p, 1);
    p /= psum;
    L = _GlintLevels4;
}

void GetPreIntegratedLightLevelsFromSky(float3 iblR, float iblPerceptualRoughness, out float2x4 L, out float2x4 p, int sliceIndex = 0)
{
    float iblMipLevel = PerceptualRoughnessToMipmapLevel(iblPerceptualRoughness);
    p[0] = SAMPLE_TEXTURECUBE_ARRAY_LOD(_SkyTextureGlintLevelWeights, s_trilinear_clamp_sampler, iblR, sliceIndex*GLINT_LEVELS_SLICE_COUNT + 0, iblMipLevel);
    p[1] = SAMPLE_TEXTURECUBE_ARRAY_LOD(_SkyTextureGlintLevelWeights, s_trilinear_clamp_sampler, iblR, sliceIndex*GLINT_LEVELS_SLICE_COUNT + 1, iblMipLevel);
    // NOTE: Normalization is required for low-precision texture formats only.
    float psum = dot(p[0], 1) + dot(p[1], 1);
    p /= psum;
    L = float2x4(_GlintLevels8_0, _GlintLevels8_1);
}


#ifdef HAS_LIGHTLOOP
void GetPreIntegratedLightLevelsFromEnv(LightLoopContext lightLoopContext, float3 iblR, float iblPerceptualRoughness, out float4 L, out float4 p, int sliceIdx = 0)
{
    // C.f. float4 SampleEnv(LightLoopContext lightLoopContext, int index, float3 texCoord, float lod, float rangeCompressionFactorCompensation, float2 positionNDC, int sliceIdx = 0) in LightLoopDef.hlsl
    if (lightLoopContext.sampleReflection != SINGLE_PASS_CONTEXT_SAMPLE_SKY)
    {
        // We only support sky reflections for now!
        L = 0.0.rrrr;
        p = 0.0.rrrr;
        return;
    }
    GetPreIntegratedLightLevelsFromSky(iblR, iblPerceptualRoughness, L, p, sliceIdx);
}
void GetPreIntegratedLightLevelsFromEnv(LightLoopContext lightLoopContext, float3 iblR, float iblPerceptualRoughness, out float2x4 L, out float2x4 p, int sliceIdx = 0)
{
    // C.f. float4 SampleEnv(LightLoopContext lightLoopContext, int index, float3 texCoord, float lod, float rangeCompressionFactorCompensation, float2 positionNDC, int sliceIdx = 0) in LightLoopDef.hlsl
    if (lightLoopContext.sampleReflection != SINGLE_PASS_CONTEXT_SAMPLE_SKY)
    {
        // We only support sky reflections for now!
        L = float2x4(0.0.rrrr, 0.0.rrrr);
        p = float2x4(0.0.rrrr, 0.0.rrrr);
        return;
    }
    GetPreIntegratedLightLevelsFromSky(iblR, iblPerceptualRoughness, L, p, sliceIdx);
}
#endif // HAS_LIGHTLOOP


#endif // PREINTEGRATED_LIGHT_DIR_HLSL
