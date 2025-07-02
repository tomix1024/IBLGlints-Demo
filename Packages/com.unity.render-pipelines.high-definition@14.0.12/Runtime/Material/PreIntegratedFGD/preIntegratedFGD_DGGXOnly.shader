Shader "Hidden/HDRP/preIntegratedFGD_DGGXOnly"
{
    SubShader
    {
        Tags{ "RenderPipeline" = "HDRenderPipeline" }
        Pass
        {
            ZTest Always Cull Off ZWrite Off

            HLSLPROGRAM

            #pragma editor_sync_compilation

            #pragma vertex Vert
            #pragma fragment Frag
            #pragma target 4.5
            #pragma only_renderers d3d11 playstation xboxone xboxseries vulkan metal switch
            #define PREFER_HALF 0
            #include "Packages/com.unity.render-pipelines.core/ShaderLibrary/Common.hlsl"
            #include "Packages/com.unity.render-pipelines.core/ShaderLibrary/ImageBasedLighting.hlsl"
            #include "Packages/com.unity.render-pipelines.high-definition/Runtime/ShaderLibrary/ShaderVariables.hlsl"
            #include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/PreIntegratedFGD/PreIntegratedFGD.cs.hlsl"



            real IntegrateDGGXOnlyFGD(real NdotV, real roughness, uint sampleCount = 4096)
            {
                // Note that our LUT covers the full [0, 1] range.
                // Therefore, we don't really want to clamp NdotV here (else the lerp slope is wrong).
                // However, if NdotV is 0, the integral is 0, so that's not what we want, either.
                // Our runtime NdotV bias is quite large, so we use a smaller one here instead.
                NdotV     = max(NdotV, REAL_EPS);
                real3 V   = real3(sqrt(1 - NdotV * NdotV), 0, NdotV);
                real acc  = 0.0;

                real3x3 localToWorld = k_identity3x3;

                for (uint i = 0; i < sampleCount; ++i)
                {
                    real2 u = Hammersley2d(i, sampleCount);

                    real3 L; // unused
                    real NdotL;
                    real NdotH;
                    real VdotH; // unused
                    SampleGGXDir(u, V, localToWorld, roughness, L, NdotL, NdotH, VdotH);
                    // pdf(H) = D_GGX * NdotH
                    // pdf(L) = D_GGX * NdotH / (4*VdotH)
                    // integrand(H) = D_GGX
                    // integrand(L) = D_GGX / (4*VdotH)

                    if (NdotL > 0.0)
                    {
                        real weightOverPdf = rcp(NdotH);
                        acc += weightOverPdf;
                    }
                }

                acc /= sampleCount;
                return acc;
            }



            struct Attributes
            {
                uint vertexID : SV_VertexID;
            };

            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float2 texCoord   : TEXCOORD0;
            };

            Varyings Vert(Attributes input)
            {
                Varyings output;

                output.positionCS = GetFullScreenTriangleVertexPosition(input.vertexID);
                output.texCoord   = GetFullScreenTriangleTexCoord(input.vertexID);

                return output;
            }

            float4 Frag(Varyings input) : SV_Target
            {
                // We want the LUT to contain the entire [0, 1] range, without losing half a texel at each side.
                float2 coordLUT = RemapHalfTexelCoordTo01(input.texCoord, FGDTEXTURE_RESOLUTION);

                // The FGD texture is parametrized as follows:
                // X = sqrt(dot(N, V))
                // Y = perceptualRoughness
                // These coordinate sampling must match the decoding in GetPreIntegratedDFG in Lit.hlsl,
                // i.e here we use perceptualRoughness, must be the same in shader
                // Note: with this angular parametrization, the LUT is almost perfectly linear,
                // except for the grazing angle when (NdotV -> 0).
                float NdotV = coordLUT.x * coordLUT.x;
                float perceptualRoughness = coordLUT.y;

                // Pre integrate GGX with smithJoint visibility as well as DisneyDiffuse
                float preFGD = IntegrateDGGXOnlyFGD(NdotV, PerceptualRoughnessToRoughness(perceptualRoughness));

                // Texture stores values in [0, 1] but represents values in [0, 2].
                return float4(0.5*preFGD, 0, 0, 1.0);
            }

            ENDHLSL
        }
    }
    Fallback Off
}
