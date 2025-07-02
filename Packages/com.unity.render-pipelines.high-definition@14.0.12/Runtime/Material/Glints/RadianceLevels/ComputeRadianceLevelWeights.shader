Shader "Hidden/HDRP/ComputeRadianceLevelWeights"
{
    SubShader
    {
        Tags{ "RenderPipeline" = "HDRenderPipeline" }
        Pass
        {
            Cull   Off
            ZTest  Off
            ZWrite Off
            Blend  Off

            HLSLPROGRAM
            #pragma editor_sync_compilation
            #pragma target 4.5
            #pragma only_renderers d3d11 playstation xboxone xboxseries vulkan metal switch

            #pragma vertex Vert
            #pragma fragment Frag

            #include "Packages/com.unity.render-pipelines.core/ShaderLibrary/Common.hlsl"
            #include "Packages/com.unity.render-pipelines.high-definition/Runtime/ShaderLibrary/ShaderVariables.hlsl"

            TEXTURECUBE(_MainTex);
            SAMPLER(sampler_MainTex);

            float4x4 _PixelCoordToViewDirWS; // Actually just 3x3, but Unity can only set 4x4

            int _OutputIndex; // This is a workaround for not fully supported MRT...

            int _SmoothStepCount; // How often to apply the smooth step function?
            int _GlintLogDistribution; // Should glints be distributed in log space across the levels?

            cbuffer _GlintLevelsData{
                float4 _GlintLevels8_0;
                float4 _GlintLevels8_1;
                float4 _GlintLevels4;
                float4 _GlintLevelsDummy;
            };

            struct Attributes
            {
                uint vertexID : SV_VertexID;
            };

            struct Varyings
            {
                float4 positionCS : SV_POSITION;
            };


            float3 S1(float3 x)
            {
                x = (-2*x + 3)*x*x;
                x = clamp(x, 0, 1);
                return x;
            }


            Varyings Vert(Attributes input)
            {
                Varyings output;
                output.positionCS = GetFullScreenTriangleVertexPosition(input.vertexID);
                return output;
            }

            float4 Frag(Varyings input) : SV_Target0
                /*out float4 output0 : SV_Target0,
                out float4 output1 : SV_Target1,
                out float4 output2 : SV_Target2)*/
            {
                // Points towards the camera
                float3 dirWS = -normalize(mul(float3(input.positionCS.xy, 1.0), (float3x3)_PixelCoordToViewDirWS));
                // Reverse it to point into the scene

                // NO LOD! want to filter across (larger) texture (?) // s_linear_clamp_sampler
                real3 val = SAMPLE_TEXTURECUBE(_MainTex, sampler_MainTex, dirWS).rgb;

                float weight = dot(val, real3(0.2126, 0.7152, 0.0722));

                float4 levels4 = _GlintLevels4;

                float3 alpha4 = ((weight - levels4.xyz) / max(levels4.yzw - levels4.xyz, 1e-12));
                bool3 mask4 = (weight <= levels4.yzw);
                mask4.yz = mask4.yz * (weight > levels4.yz);

                if (_GlintLogDistribution)
                {
                    // We know that level0 = 0.
                    float4 logLevels4 = log(levels4);
                    float logWeight = log(max(levels4.y, weight));
                    float3 logAlpha4 = ((logWeight - logLevels4.xyz) / (logLevels4.yzw - logLevels4.xyz));
                    alpha4.yz = logAlpha4.yz; // First interpolation is linear here!
                }
                for (int i = 0; i < min(_SmoothStepCount, 2); ++i)
                    alpha4 = S1(alpha4);

                float4 discreteValues4 = float4(
                    (1-alpha4.x)*mask4.x,
                    (1-alpha4.yz)*mask4.yz + alpha4.xy * mask4.xy,
                    alpha4.z * mask4.z
                );

                float2x4 levels8 = float2x4(_GlintLevels8_0, _GlintLevels8_1);

                float3 alpha8_0 = ((weight - levels8[0].xyz) / max(levels8[0].yzw - levels8[0].xyz, 1e-12));
                float alpha8_01 = ((weight - levels8[0].w)   / max(levels8[1].x   - levels8[0].w,   1e-12));
                float3 alpha8_1 = ((weight - levels8[1].xyz) / max(levels8[1].yzw - levels8[1].xyz, 1e-12));

                bool3 mask8_0 = (weight <= levels8[0].yzw);
                mask8_0.yz = mask8_0.yz * (weight > levels8[0].yz);
                bool mask8_01 = (weight > levels8[0].w) * (weight <= levels8[1].x);
                bool3 mask8_1 = (weight > levels8[1].xyz) * (weight <= levels8[1].yzw);

                if (_GlintLogDistribution)
                {
                    // We know that level0 = 0.
                    float2x4 logLevels8 = log(levels8);
                    float logWeight = log(max(levels8[0].y, weight));
                    float3 logAlpha8_0 = ((logWeight - logLevels8[0].xyz) / (logLevels8[0].yzw - logLevels8[0].xyz));
                    float logAlpha8_01 = ((logWeight - logLevels8[0].w) / (logLevels8[1].x - logLevels8[0].w));
                    float3 logAlpha8_1 = ((logWeight - logLevels8[1].xyz) / (logLevels8[1].yzw - logLevels8[1].xyz));
                    alpha8_0.yz = logAlpha8_0.yz; // First interpolation is linear here!
                    alpha8_01 = logAlpha8_01;
                    alpha8_1 = logAlpha8_1;
                }
                for (int j = 0; j < min(_SmoothStepCount, 2); ++j)
                {
                    alpha8_0 = S1(alpha8_0);
                    alpha8_01 = S1(alpha8_01.xxx).x;
                    alpha8_1 = S1(alpha8_1);
                }

                float2x4 discreteValues8 = float2x4(
                    float4(
                        (1-alpha8_0.x)*mask8_0.x,
                        (1-alpha8_0.yz)*mask8_0.yz + alpha8_0.xy * mask8_0.xy,
                        (1-alpha8_01)*mask8_01 + alpha8_0.z * mask8_0.z
                    ),
                    float4(
                        (1-alpha8_1.x)*mask8_1.x + alpha8_01 * mask8_01,
                        (1-alpha8_1.yz)*mask8_1.yz + alpha8_1.xy * mask8_1.xy,
                        alpha8_1.z * mask8_1.z
                    )
                );

                switch (_OutputIndex)
                {
                    case 0:
                        return discreteValues8[0];
                    case 1:
                        return discreteValues8[1];
                    case 2:
                        return discreteValues4;
                    default:
                        return float4(1,0,1,1);
                }

                // output0 = float4(1,0,1,1); // float4(X_Y_Z, weight);
                // output1 = float4(1,0,1,1); // float4(X2_Y2_Z2, weight);
                // output2 = float4(1,0,1,1); // float4(XY_YZ_ZX, weight);
            }
            ENDHLSL
        }
    }
}
