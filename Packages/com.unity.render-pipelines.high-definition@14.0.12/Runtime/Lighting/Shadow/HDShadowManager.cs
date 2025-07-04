using System;
using System.Collections.Generic;
using UnityEngine.Experimental.Rendering;
using UnityEngine.Experimental.Rendering.RenderGraphModule;
using UnityEngine.Serialization;

namespace UnityEngine.Rendering.HighDefinition
{
    /// <summary>
    /// Shadow Filtering Quality
    /// </summary>
    public enum HDShadowFilteringQuality
    {
        /// <summary>
        /// Low Shadow Filtering Quality
        /// </summary>
        Low = 0,
        /// <summary>
        /// Medium Shadow Filtering Quality
        /// </summary>
        Medium = 1,
        /// <summary>
        /// High Shadow Filtering Quality
        /// </summary>
        High = 2
    }

    /// <summary>
    /// Area Shadow Filtering Quality
    /// </summary>
    public enum HDAreaShadowFilteringQuality
    {
        /// <summary>
        /// Area Medium Shadow Filtering Quality
        /// </summary>
        Medium = 0,
        /// <summary>
        /// Area High Shadow Filtering Quality
        /// </summary>
        High = 1
    }

    enum ShadowMapType
    {
        CascadedDirectional,
        PunctualAtlas,
        AreaLightAtlas
    }

    enum ShadowMapUpdateType
    {
        // Fully dynamic shadow maps
        Dynamic = 0,
        // Fully cached shadow maps (nothing is rendered unless requested)
        Cached,
        // Mixed, static shadow caster are cached and updated as indicated, dynamic are drawn on top.
        Mixed
    }

    [GenerateHLSL(needAccessors = false)]
    struct HDShadowData
    {
        public Vector3 rot0;
        public Vector3 rot1;
        public Vector3 rot2;
        public Vector3 pos;
        public Vector4 proj;

        public Vector2 atlasOffset;
        public float worldTexelSize;
        public float normalBias;

        [SurfaceDataAttributes(precision = FieldPrecision.Real)]
        public Vector4 zBufferParam;
        public Vector4 shadowMapSize;

        public Vector4 shadowFilterParams0;

        public Vector3 cacheTranslationDelta;
        public float isInCachedAtlas;

        public Matrix4x4 shadowToWorld;
    }

    // We use a different structure for directional light because these is a lot of data there
    // and it will add too much useless stuff for other lights
    // Note: In order to support HLSL array generation, we need to use fixed arrays and so a unsafe context for this struct
    [GenerateHLSL(needAccessors = false)]
    unsafe struct HDDirectionalShadowData
    {
        // We can't use Vector4 here because the vector4[] makes this struct non blittable
        [HLSLArray(4, typeof(Vector4))]
        public fixed float sphereCascades[4 * 4];

        [SurfaceDataAttributes(precision = FieldPrecision.Real)]
        public Vector4 cascadeDirection;

        [HLSLArray(4, typeof(float))]
        [SurfaceDataAttributes(precision = FieldPrecision.Real)]
        public fixed float cascadeBorders[4];

        public float fadeScale;
        public float fadeBias;
    }

    class HDShadowRequest
    {
        public Matrix4x4 view;
        // Use the y flipped device projection matrix as light projection matrix
        public Matrix4x4 deviceProjectionYFlip;
        public Matrix4x4 deviceProjection;
        public Matrix4x4 projection;
        public BatchCullingProjectionType projectionType;
        public Matrix4x4 shadowToWorld;
        public Vector3 position;
        public Vector4 zBufferParam;
        // Warning: these viewport fields are updated by ProcessShadowRequests and are invalid before
        public Rect dynamicAtlasViewport;
        public Rect cachedAtlasViewport;
        public bool zClip;
        public Vector4[] frustumPlanes;

        // Store the final shadow indice in the shadow data array
        // Warning: the index is computed during ProcessShadowRequest and so is invalid before calling this function
        public int shadowIndex;

        // Determine in which atlas the shadow will be rendered
        public ShadowMapType shadowMapType = ShadowMapType.PunctualAtlas;

        // TODO: Remove these field once scriptable culling is here (currently required by ScriptableRenderContext.DrawShadows)
        public int lightIndex;
        public ShadowSplitData splitData;
        // end

        public float normalBias;
        public float worldTexelSize;
        public float slopeBias;

        // PCSS parameters
        public float shadowSoftness;
        public int blockerSampleCount;
        public int filterSampleCount;
        public float minFilterSize;

        // IMS parameters
        public float kernelSize;
        public float lightAngle;
        public float maxDepthBias;

        public Vector4 evsmParams;

        public bool shouldUseCachedShadowData = false;
        public bool shouldRenderCachedComponent = false;

        public HDShadowData cachedShadowData;

        public bool isInCachedAtlas;
        public bool isMixedCached = false;
    }

    enum DirectionalShadowAlgorithm
    {
        PCF5x5,
        PCF7x7,
        PCSS,
        IMS
    }

    /// <summary>
    /// Screen Space Shadows format.
    /// </summary>
    public enum ScreenSpaceShadowFormat
    {
        /// <summary>R8G8B8A8 format for fastest rendering.</summary>
        R8G8B8A8 = GraphicsFormat.R8G8B8A8_UNorm,
        /// <summary>R16G16B16A16 format for better quality.</summary>
        R16G16B16A16 = GraphicsFormat.R16G16B16A16_SFloat
    }

    /// <summary>
    /// Shadows Global Settings.
    /// </summary>
    [Serializable]
    public struct HDShadowInitParameters
    {
        /// <summary>
        /// Shadow Atlases parameters.
        /// </summary>
        [Serializable]
        public struct HDShadowAtlasInitParams
        {
            /// <summary>Shadow Atlas resolution.</summary>
            public int shadowAtlasResolution;
            /// <summary>Shadow Atlas depth bits.</summary>
            public DepthBits shadowAtlasDepthBits;
            /// <summary>Enable dynamic rescale of the atlas.</summary>
            public bool useDynamicViewportRescale;

            internal static HDShadowAtlasInitParams GetDefault()
            {
                return new HDShadowAtlasInitParams()
                {
                    shadowAtlasResolution = k_DefaultShadowAtlasResolution,
                    shadowAtlasDepthBits = k_DefaultShadowMapDepthBits,
                    useDynamicViewportRescale = true
                };
            }
        }

        internal static HDShadowInitParameters NewDefault() => new HDShadowInitParameters()
        {
            maxShadowRequests = k_DefaultMaxShadowRequests,
            directionalShadowsDepthBits = k_DefaultShadowMapDepthBits,
            punctualLightShadowAtlas = HDShadowAtlasInitParams.GetDefault(),
            areaLightShadowAtlas = HDShadowAtlasInitParams.GetDefault(),
            cachedPunctualLightShadowAtlas = 2048,
            cachedAreaLightShadowAtlas = 1024,
            allowDirectionalMixedCachedShadows = false,
            shadowResolutionDirectional = new IntScalableSetting(new[] { 256, 512, 1024, 2048 }, ScalableSettingSchemaId.With4Levels),
            shadowResolutionArea = new IntScalableSetting(new[] { 256, 512, 1024, 2048 }, ScalableSettingSchemaId.With4Levels),
            shadowResolutionPunctual = new IntScalableSetting(new[] { 256, 512, 1024, 2048 }, ScalableSettingSchemaId.With4Levels),
            shadowFilteringQuality = HDShadowFilteringQuality.Medium,
            areaShadowFilteringQuality = HDAreaShadowFilteringQuality.Medium,
            supportScreenSpaceShadows = false,
            maxScreenSpaceShadowSlots = 4,
            screenSpaceShadowBufferFormat = ScreenSpaceShadowFormat.R16G16B16A16,
            maxDirectionalShadowMapResolution = 2048,
            maxAreaShadowMapResolution = 2048,
            maxPunctualShadowMapResolution = 2048,
        };

        internal const int k_DefaultShadowAtlasResolution = 4096;
        internal const int k_DefaultMaxShadowRequests = 128;
        internal const DepthBits k_DefaultShadowMapDepthBits = DepthBits.Depth32;

        /// <summary>Maximum number of shadow requests at the same time.</summary>
        public int maxShadowRequests;
        /// <summary>Depth bits for directional shadows.</summary>
        public DepthBits directionalShadowsDepthBits;

        /// <summary>Shadow filtering quality.</summary>
        [FormerlySerializedAs("shadowQuality")]
        public HDShadowFilteringQuality shadowFilteringQuality;

        /// <summary>Area Shadow filtering quality.</summary>
        public HDAreaShadowFilteringQuality areaShadowFilteringQuality;

        /// <summary>Initialization parameters for punctual shadows atlas.</summary>
        public HDShadowAtlasInitParams punctualLightShadowAtlas;
        /// <summary>Initialization parameters for area shadows atlas.</summary>
        public HDShadowAtlasInitParams areaLightShadowAtlas;

        /// <summary>Resolution for the punctual lights cached shadow maps atlas.</summary>
        public int cachedPunctualLightShadowAtlas;

        /// <summary>Resolution for the area lights cached shadow maps atlas.</summary>
        public int cachedAreaLightShadowAtlas;

        /// <summary>Maximum shadow map resolution for directional lights.</summary>
        public bool allowDirectionalMixedCachedShadows;


        /// <summary>Shadow scalable resolution for directional lights.</summary>
        public IntScalableSetting shadowResolutionDirectional;
        /// <summary>Shadow scalable resolution for point lights.</summary>
        public IntScalableSetting shadowResolutionPunctual;
        /// <summary>Shadow scalable resolution for area lights.</summary>
        public IntScalableSetting shadowResolutionArea;

        /// <summary>Maximum shadow map resolution for directional lights.</summary>
        public int maxDirectionalShadowMapResolution;
        /// <summary>Maximum shadow map resolution for punctual lights.</summary>
        public int maxPunctualShadowMapResolution;
        /// <summary>Maximum shadow map resolution for area lights.</summary>
        public int maxAreaShadowMapResolution;

        /// <summary>Enable support for screen space shadows.</summary>
        public bool supportScreenSpaceShadows;
        /// <summary>Maximum number of screen space shadows.</summary>
        public int maxScreenSpaceShadowSlots;
        /// <summary>Format for screen space shadows.</summary>
        public ScreenSpaceShadowFormat screenSpaceShadowBufferFormat;
    }

    class HDShadowResolutionRequest
    {
        public Rect dynamicAtlasViewport;
        public Rect cachedAtlasViewport;
        public Vector2 resolution;
        public ShadowMapType shadowMapType;

        public HDShadowResolutionRequest ShallowCopy()
        {
            return (HDShadowResolutionRequest)this.MemberwiseClone();
        }
    }

    internal struct ShadowResult
    {
        public TextureHandle punctualShadowResult;
        public TextureHandle cachedPunctualShadowResult;
        public TextureHandle directionalShadowResult;
        public TextureHandle areaShadowResult;
        public TextureHandle cachedAreaShadowResult;
    }


    class HDShadowManager
    {
        public const int k_DirectionalShadowCascadeCount = 4;
        public const int k_MinShadowMapResolution = 16;
        public const int k_MaxShadowMapResolution = 16384;

        List<HDShadowData> m_ShadowDatas = new List<HDShadowData>();
        HDShadowRequest[] m_ShadowRequests;
        HDShadowResolutionRequest[] m_ShadowResolutionRequests;
        HDDirectionalShadowData[] m_CachedDirectionalShadowData;

        HDDirectionalShadowData m_DirectionalShadowData;

        // Structured buffer of shadow datas
        ComputeBuffer m_ShadowDataBuffer;
        ComputeBuffer m_DirectionalShadowDataBuffer;

        // The two shadowmaps atlases we uses, one for directional cascade (without resize) and the second for the rest of the shadows
        HDDynamicShadowAtlas m_CascadeAtlas;
        HDDynamicShadowAtlas m_Atlas;
        HDDynamicShadowAtlas m_AreaLightShadowAtlas;

        int m_MaxShadowRequests;
        int m_ShadowRequestCount;
        int m_CascadeCount;
        int m_ShadowResolutionRequestCounter;

        Material m_ClearShadowMaterial;
        Material m_BlitShadowMaterial;

        ConstantBuffer<ShaderVariablesGlobal> m_GlobalShaderVariables;

        public static HDCachedShadowManager cachedShadowManager { get { return HDCachedShadowManager.instance; } }

        public void InitShadowManager(HDRenderPipelineRuntimeResources renderPipelineResources, HDShadowInitParameters initParams, RenderGraph renderGraph, Shader clearShader)
        {
            // Even when shadows are disabled (maxShadowRequests == 0) we need to allocate compute buffers to avoid having
            // resource not bound errors when dispatching a compute shader.
            m_ShadowDataBuffer = new ComputeBuffer(Mathf.Max(initParams.maxShadowRequests, 1), System.Runtime.InteropServices.Marshal.SizeOf(typeof(HDShadowData)));
            m_DirectionalShadowDataBuffer = new ComputeBuffer(1, System.Runtime.InteropServices.Marshal.SizeOf(typeof(HDDirectionalShadowData)));
            m_MaxShadowRequests = initParams.maxShadowRequests;
            m_ShadowRequestCount = 0;

            if (initParams.maxShadowRequests == 0)
                return;

            m_ClearShadowMaterial = CoreUtils.CreateEngineMaterial(clearShader);
            m_BlitShadowMaterial = CoreUtils.CreateEngineMaterial(renderPipelineResources.shaders.shadowBlitPS);

            // Prevent the list from resizing their internal container when we add shadow requests
            m_ShadowDatas.Capacity = Math.Max(initParams.maxShadowRequests, m_ShadowDatas.Capacity);
            m_ShadowResolutionRequests = new HDShadowResolutionRequest[initParams.maxShadowRequests];
            m_ShadowRequests = new HDShadowRequest[initParams.maxShadowRequests];
            m_CachedDirectionalShadowData = new HDDirectionalShadowData[1]; // we only support directional light shadow

            m_GlobalShaderVariables = new ConstantBuffer<ShaderVariablesGlobal>();

            for (int i = 0; i < initParams.maxShadowRequests; i++)
            {
                m_ShadowResolutionRequests[i] = new HDShadowResolutionRequest();
            }

            HDShadowAtlas.HDShadowAtlasInitParameters punctualAtlasInitParams = new HDShadowAtlas.HDShadowAtlasInitParameters(renderPipelineResources, renderGraph, useSharedTexture: false, initParams.punctualLightShadowAtlas.shadowAtlasResolution,
                initParams.punctualLightShadowAtlas.shadowAtlasResolution, HDShaderIDs._ShadowmapAtlas, m_ClearShadowMaterial, initParams.maxShadowRequests, initParams, m_GlobalShaderVariables);
            punctualAtlasInitParams.name = "Shadow Map Atlas";

            // The cascade atlas will be allocated only if there is a directional light
            m_Atlas = new HDDynamicShadowAtlas(punctualAtlasInitParams);
            // Cascade atlas render texture will only be allocated if there is a shadow casting directional light
            HDShadowAtlas.BlurAlgorithm cascadeBlur = GetDirectionalShadowAlgorithm() == DirectionalShadowAlgorithm.IMS ? HDShadowAtlas.BlurAlgorithm.IM : HDShadowAtlas.BlurAlgorithm.None;

            HDShadowAtlas.HDShadowAtlasInitParameters dirAtlasInitParams = punctualAtlasInitParams;
            dirAtlasInitParams.useSharedTexture = true;
            dirAtlasInitParams.width = 1;
            dirAtlasInitParams.height = 1;
            dirAtlasInitParams.atlasShaderID = HDShaderIDs._ShadowmapCascadeAtlas;
            dirAtlasInitParams.blurAlgorithm = cascadeBlur;
            dirAtlasInitParams.depthBufferBits = initParams.directionalShadowsDepthBits;
            dirAtlasInitParams.name = "Cascade Shadow Map Atlas";

            m_CascadeAtlas = new HDDynamicShadowAtlas(dirAtlasInitParams);

            HDShadowAtlas.HDShadowAtlasInitParameters areaAtlasInitParams = punctualAtlasInitParams;
            if (ShaderConfig.s_AreaLights == 1)
            {
                areaAtlasInitParams.useSharedTexture = false;
                areaAtlasInitParams.width = initParams.areaLightShadowAtlas.shadowAtlasResolution;
                areaAtlasInitParams.height = initParams.areaLightShadowAtlas.shadowAtlasResolution;
                areaAtlasInitParams.atlasShaderID = HDShaderIDs._ShadowmapAreaAtlas;
                areaAtlasInitParams.blurAlgorithm = GetAreaLightShadowBlurAlgorithm();
                areaAtlasInitParams.depthBufferBits = initParams.areaLightShadowAtlas.shadowAtlasDepthBits;
                areaAtlasInitParams.name = "Area Light Shadow Map Atlas";


                m_AreaLightShadowAtlas = new HDDynamicShadowAtlas(areaAtlasInitParams);
            }

            HDShadowAtlas.HDShadowAtlasInitParameters cachedPunctualAtlasInitParams = punctualAtlasInitParams;
            cachedPunctualAtlasInitParams.useSharedTexture = true;
            cachedPunctualAtlasInitParams.width = initParams.cachedPunctualLightShadowAtlas;
            cachedPunctualAtlasInitParams.height = initParams.cachedPunctualLightShadowAtlas;
            cachedPunctualAtlasInitParams.atlasShaderID = HDShaderIDs._CachedShadowmapAtlas;
            cachedPunctualAtlasInitParams.name = "Cached Shadow Map Atlas";
            cachedPunctualAtlasInitParams.isShadowCache = true;

            cachedShadowManager.InitPunctualShadowAtlas(cachedPunctualAtlasInitParams);
            if (ShaderConfig.s_AreaLights == 1)
            {
                HDShadowAtlas.HDShadowAtlasInitParameters cachedAreaAtlasInitParams = areaAtlasInitParams;
                cachedAreaAtlasInitParams.useSharedTexture = true;
                cachedAreaAtlasInitParams.width = initParams.cachedAreaLightShadowAtlas;
                cachedAreaAtlasInitParams.height = initParams.cachedAreaLightShadowAtlas;
                cachedAreaAtlasInitParams.atlasShaderID = HDShaderIDs._CachedAreaLightShadowmapAtlas;
                cachedAreaAtlasInitParams.name = "Cached Area Light Shadow Map Atlas";
                cachedAreaAtlasInitParams.isShadowCache = true;

                cachedShadowManager.InitAreaLightShadowAtlas(cachedAreaAtlasInitParams);
            }

            cachedShadowManager.InitDirectionalState(dirAtlasInitParams, initParams.allowDirectionalMixedCachedShadows);
        }

        public void Cleanup(RenderGraph renderGraph)
        {
            m_ShadowDataBuffer.Dispose();
            m_DirectionalShadowDataBuffer.Dispose();

            if (m_MaxShadowRequests == 0)
                return;

            m_Atlas.Release(renderGraph);
            if (ShaderConfig.s_AreaLights == 1)
                m_AreaLightShadowAtlas.Release(renderGraph);
            m_CascadeAtlas.Release(renderGraph);

            CoreUtils.Destroy(m_ClearShadowMaterial);
            cachedShadowManager.Cleanup(renderGraph);

            m_GlobalShaderVariables.Release();
        }

        // Keep in sync with both HDShadowSampling.hlsl
        public static DirectionalShadowAlgorithm GetDirectionalShadowAlgorithm()
        {
            switch (HDRenderPipeline.currentAsset.currentPlatformRenderPipelineSettings.hdShadowInitParams.shadowFilteringQuality)
            {
                case HDShadowFilteringQuality.Low:
                {
                    return DirectionalShadowAlgorithm.PCF5x5;
                }
                case HDShadowFilteringQuality.Medium:
                {
                    return DirectionalShadowAlgorithm.PCF7x7;
                }
                case HDShadowFilteringQuality.High:
                {
                    return DirectionalShadowAlgorithm.PCSS;
                }
            }
            ;
            return DirectionalShadowAlgorithm.PCF5x5;
        }

        public static HDShadowAtlas.BlurAlgorithm GetAreaLightShadowBlurAlgorithm()
        {
            return HDRenderPipeline.currentAsset.currentPlatformRenderPipelineSettings.hdShadowInitParams.areaShadowFilteringQuality == HDAreaShadowFilteringQuality.High ?
                HDShadowAtlas.BlurAlgorithm.None : HDShadowAtlas.BlurAlgorithm.EVSM;
        }

        public void UpdateShaderVariablesGlobalCB(ref ShaderVariablesGlobal cb)
        {
            if (m_MaxShadowRequests == 0)
                return;

            cb._CascadeShadowCount = (uint)(m_CascadeCount + 1);
            cb._ShadowAtlasSize = new Vector4(m_Atlas.width, m_Atlas.height, 1.0f / m_Atlas.width, 1.0f / m_Atlas.height);
            cb._CascadeShadowAtlasSize = new Vector4(m_CascadeAtlas.width, m_CascadeAtlas.height, 1.0f / m_CascadeAtlas.width, 1.0f / m_CascadeAtlas.height);
            cb._CachedShadowAtlasSize = new Vector4(cachedShadowManager.punctualShadowAtlas.width, cachedShadowManager.punctualShadowAtlas.height, 1.0f / cachedShadowManager.punctualShadowAtlas.width, 1.0f / cachedShadowManager.punctualShadowAtlas.height);
            if (ShaderConfig.s_AreaLights == 1)
            {
                cb._AreaShadowAtlasSize = new Vector4(m_AreaLightShadowAtlas.width, m_AreaLightShadowAtlas.height, 1.0f / m_AreaLightShadowAtlas.width, 1.0f / m_AreaLightShadowAtlas.height);
                cb._CachedAreaShadowAtlasSize = new Vector4(cachedShadowManager.areaShadowAtlas.width, cachedShadowManager.areaShadowAtlas.height, 1.0f / cachedShadowManager.areaShadowAtlas.width, 1.0f / cachedShadowManager.areaShadowAtlas.height);
            }
        }

        public void UpdateDirectionalShadowResolution(int resolution, int cascadeCount)
        {
            Vector2Int atlasResolution = new Vector2Int(resolution, resolution);

            if (cascadeCount > 1)
                atlasResolution.x *= 2;
            if (cascadeCount > 2)
                atlasResolution.y *= 2;

            m_CascadeAtlas.UpdateSize(atlasResolution);
            if (cachedShadowManager.DirectionalHasCachedAtlas())
                cachedShadowManager.directionalLightAtlas.UpdateSize(atlasResolution);
        }

        internal int ReserveShadowResolutions(Vector2 resolution, ShadowMapType shadowMapType, int lightID, int index, ShadowMapUpdateType updateType)
        {
            if (m_ShadowRequestCount >= m_MaxShadowRequests)
            {
#if UNITY_EDITOR || DEVELOPMENT_BUILD
                Debug.LogWarning("Max shadow requests count reached, dropping all exceeding requests. You can increase this limit by changing the Maximum Shadows on Screen property in the HDRP asset.");
#endif
                return -1;
            }

            m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter].shadowMapType = shadowMapType;

            // Note: for cached shadows we manage the resolution requests directly on the CachedShadowAtlas as they need special handling. We however keep incrementing the counter for two reasons:
            //      - Maintain the limit of m_MaxShadowRequests
            //      - Avoid to refactor other parts that the shadow manager that get requests indices from here.

            if (updateType != ShadowMapUpdateType.Cached || shadowMapType == ShadowMapType.CascadedDirectional)
            {
                m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter].resolution = resolution;
                m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter].dynamicAtlasViewport.width = resolution.x;
                m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter].dynamicAtlasViewport.height = resolution.y;

                switch (shadowMapType)
                {
                    case ShadowMapType.PunctualAtlas:
                        m_Atlas.ReserveResolution(m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter]);
                        break;
                    case ShadowMapType.AreaLightAtlas:
                        m_AreaLightShadowAtlas.ReserveResolution(m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter]);
                        break;
                    case ShadowMapType.CascadedDirectional:
                        m_CascadeAtlas.ReserveResolution(m_ShadowResolutionRequests[m_ShadowResolutionRequestCounter]);
                        break;
                }
            }


            m_ShadowResolutionRequestCounter++;
            m_ShadowRequestCount = m_ShadowResolutionRequestCounter;

            return m_ShadowResolutionRequestCounter - 1;
        }

        internal HDShadowResolutionRequest GetResolutionRequest(int index)
        {
            if (index < 0 || index >= m_ShadowRequestCount)
                return null;

            return m_ShadowResolutionRequests[index];
        }

        public Vector2 GetReservedResolution(int index)
        {
            if (index < 0 || index >= m_ShadowRequestCount)
                return Vector2.zero;

            return m_ShadowResolutionRequests[index].resolution;
        }

        internal void UpdateShadowRequest(int index, HDShadowRequest shadowRequest, ShadowMapUpdateType updateType)
        {
            if (index >= m_ShadowRequestCount)
                return;

            m_ShadowRequests[index] = shadowRequest;

            bool addToCached = updateType == ShadowMapUpdateType.Cached || updateType == ShadowMapUpdateType.Mixed;
            bool addDynamic = updateType == ShadowMapUpdateType.Dynamic || updateType == ShadowMapUpdateType.Mixed;

            switch (shadowRequest.shadowMapType)
            {
                case ShadowMapType.PunctualAtlas:
                {
                    if (addToCached)
                        cachedShadowManager.punctualShadowAtlas.AddShadowRequest(shadowRequest);
                    if (addDynamic)
                    {
                        m_Atlas.AddShadowRequest(shadowRequest);
                        if (updateType == ShadowMapUpdateType.Mixed)
                            m_Atlas.AddRequestToPendingBlitFromCache(shadowRequest);
                    }

                    break;
                }
                case ShadowMapType.CascadedDirectional:
                {

                    if (updateType == ShadowMapUpdateType.Mixed && cachedShadowManager.DirectionalHasCachedAtlas())    // TODO_FCC: Make sure this is disabled when setting the option off in the asset.
                    {
                        cachedShadowManager.directionalLightAtlas.AddShadowRequest(shadowRequest);
                        m_CascadeAtlas.AddRequestToPendingBlitFromCache(shadowRequest);
                    }

                    m_CascadeAtlas.AddShadowRequest(shadowRequest);
                    break;
                }
                case ShadowMapType.AreaLightAtlas:
                {
                    if (addToCached)
                    {
                        cachedShadowManager.areaShadowAtlas.AddShadowRequest(shadowRequest);
                    }
                    if (addDynamic)
                    {
                        m_AreaLightShadowAtlas.AddShadowRequest(shadowRequest);
                        if (updateType == ShadowMapUpdateType.Mixed)
                            m_AreaLightShadowAtlas.AddRequestToPendingBlitFromCache(shadowRequest);
                    }

                    break;
                }
            }
            ;
        }

        public void UpdateCascade(int cascadeIndex, Vector4 cullingSphere, float border)
        {
            if (cullingSphere.w != float.NegativeInfinity)
            {
                cullingSphere.w *= cullingSphere.w;
            }

            m_CascadeCount = Mathf.Max(m_CascadeCount, cascadeIndex);

            unsafe
            {
                fixed (float* sphereCascadesBuffer = m_DirectionalShadowData.sphereCascades)
                    ((Vector4*)sphereCascadesBuffer)[cascadeIndex] = cullingSphere;
                fixed (float* cascadeBorders = m_DirectionalShadowData.cascadeBorders)
                    cascadeBorders[cascadeIndex] = border;
            }
        }

        HDShadowData CreateShadowData(HDShadowRequest shadowRequest, HDShadowAtlas atlas)
        {
            HDShadowData data = new HDShadowData();

            var devProj = shadowRequest.deviceProjection;
            var view = shadowRequest.view;
            data.proj = new Vector4(devProj.m00, devProj.m11, devProj.m22, devProj.m23);
            data.pos = shadowRequest.position;
            data.rot0 = new Vector3(view.m00, view.m01, view.m02);
            data.rot1 = new Vector3(view.m10, view.m11, view.m12);
            data.rot2 = new Vector3(view.m20, view.m21, view.m22);
            data.shadowToWorld = shadowRequest.shadowToWorld;
            data.cacheTranslationDelta = new Vector3(0.0f, 0.0f, 0.0f);

            var viewport = shadowRequest.isInCachedAtlas ? shadowRequest.cachedAtlasViewport : shadowRequest.dynamicAtlasViewport;

            // Compute the scale and offset (between 0 and 1) for the atlas coordinates
            float rWidth = 1.0f / atlas.width;
            float rHeight = 1.0f / atlas.height;
            data.atlasOffset = Vector2.Scale(new Vector2(rWidth, rHeight), new Vector2(viewport.x, viewport.y));

            data.shadowMapSize = new Vector4(viewport.width, viewport.height, 1.0f / viewport.width, 1.0f / viewport.height);

            data.normalBias = shadowRequest.normalBias;
            data.worldTexelSize = shadowRequest.worldTexelSize;

            data.shadowFilterParams0.x = shadowRequest.shadowSoftness;
            data.shadowFilterParams0.y = HDShadowUtils.Asfloat(shadowRequest.blockerSampleCount);
            data.shadowFilterParams0.z = HDShadowUtils.Asfloat(shadowRequest.filterSampleCount);
            data.shadowFilterParams0.w = shadowRequest.minFilterSize;

            data.zBufferParam = shadowRequest.zBufferParam;
            if (atlas.HasBlurredEVSM())
            {
                data.shadowFilterParams0 = shadowRequest.evsmParams;
            }

            data.isInCachedAtlas = shadowRequest.isInCachedAtlas ? 1.0f : 0.0f;

            return data;
        }

        unsafe Vector4 GetCascadeSphereAtIndex(int index)
        {
            fixed (float* sphereCascadesBuffer = m_DirectionalShadowData.sphereCascades)
            {
                return ((Vector4*)sphereCascadesBuffer)[index];
            }
        }

        public void UpdateCullingParameters(ref ScriptableCullingParameters cullingParams, float maxShadowDistance)
        {
            cullingParams.shadowDistance = Mathf.Min(maxShadowDistance, cullingParams.shadowDistance);
        }

        public void LayoutShadowMaps(LightingDebugSettings lightingDebugSettings)
        {
            if (m_MaxShadowRequests == 0)
                return;

            cachedShadowManager.UpdateDebugSettings(lightingDebugSettings);

            m_Atlas.UpdateDebugSettings(lightingDebugSettings);

            if (m_CascadeAtlas != null)
                m_CascadeAtlas.UpdateDebugSettings(lightingDebugSettings);

            if (ShaderConfig.s_AreaLights == 1)
                m_AreaLightShadowAtlas.UpdateDebugSettings(lightingDebugSettings);

            if (lightingDebugSettings.shadowResolutionScaleFactor != 1.0f)
            {
                foreach (var shadowResolutionRequest in m_ShadowResolutionRequests)
                {
                    // We don't rescale the directional shadows with the global shadow scale factor
                    // because there is no dynamic atlas rescale when it overflow.
                    if (shadowResolutionRequest.shadowMapType != ShadowMapType.CascadedDirectional)
                        shadowResolutionRequest.resolution *= lightingDebugSettings.shadowResolutionScaleFactor;
                }
            }

            // Assign a position to all the shadows in the atlas, and scale shadows if needed
            if (m_CascadeAtlas != null && !m_CascadeAtlas.Layout(false))
                Debug.LogError("Cascade Shadow atlasing has failed, only one directional light can cast shadows at a time");
            m_Atlas.Layout();
            if (ShaderConfig.s_AreaLights == 1)
                m_AreaLightShadowAtlas.Layout();
        }

        unsafe public void PrepareGPUShadowDatas(CullingResults cullResults, HDCamera camera)
        {
            if (m_MaxShadowRequests == 0)
                return;

            int shadowIndex = 0;

            m_ShadowDatas.Clear();

            // Create all HDShadowDatas and update them with shadow request datas
            for (int i = 0; i < m_ShadowRequestCount; i++)
            {
                Debug.Assert(m_ShadowRequests[i] != null);

                HDShadowAtlas atlas = m_Atlas;
                if (m_ShadowRequests[i].isInCachedAtlas)
                {
                    atlas = cachedShadowManager.punctualShadowAtlas;
                }

                if (m_ShadowRequests[i].shadowMapType == ShadowMapType.CascadedDirectional)
                {
                    atlas = m_CascadeAtlas;
                }
                else if (m_ShadowRequests[i].shadowMapType == ShadowMapType.AreaLightAtlas)
                {
                    atlas = m_AreaLightShadowAtlas;
                    if (m_ShadowRequests[i].isInCachedAtlas)
                    {
                        atlas = cachedShadowManager.areaShadowAtlas;
                    }
                }

                HDShadowData shadowData;
                if (m_ShadowRequests[i].shouldUseCachedShadowData)
                {
                    shadowData = m_ShadowRequests[i].cachedShadowData;
                }
                else
                {
                    shadowData = CreateShadowData(m_ShadowRequests[i], atlas);
                    m_ShadowRequests[i].cachedShadowData = shadowData;
                }

                m_ShadowDatas.Add(shadowData);
                m_ShadowRequests[i].shadowIndex = shadowIndex++;
            }

            int first = k_DirectionalShadowCascadeCount, second = k_DirectionalShadowCascadeCount;

            fixed (float* sphereBuffer = m_DirectionalShadowData.sphereCascades)
            {
                Vector4* sphere = (Vector4*)sphereBuffer;
                for (int i = 0; i < k_DirectionalShadowCascadeCount; i++)
                {
                    first = (first == k_DirectionalShadowCascadeCount && sphere[i].w > 0.0f) ? i : first;
                    second = ((second == k_DirectionalShadowCascadeCount || second == first) && sphere[i].w > 0.0f) ? i : second;
                }
            }

            // Update directional datas:
            if (second != k_DirectionalShadowCascadeCount)
                m_DirectionalShadowData.cascadeDirection = (GetCascadeSphereAtIndex(second) - GetCascadeSphereAtIndex(first)).normalized;
            else
                m_DirectionalShadowData.cascadeDirection = Vector4.zero;

            HDShadowSettings shadowSettings = camera.volumeStack.GetComponent<HDShadowSettings>();
            m_DirectionalShadowData.cascadeDirection.w = shadowSettings.cascadeShadowSplitCount.value;

            GetShadowFadeScaleAndBias(shadowSettings, out m_DirectionalShadowData.fadeScale, out m_DirectionalShadowData.fadeBias);

            if (m_ShadowRequestCount > 0)
            {
                // Upload the shadow buffers to GPU
                m_ShadowDataBuffer.SetData(m_ShadowDatas);
                m_CachedDirectionalShadowData[0] = m_DirectionalShadowData;
                m_DirectionalShadowDataBuffer.SetData(m_CachedDirectionalShadowData);
            }
        }

        void GetShadowFadeScaleAndBias(HDShadowSettings shadowSettings, out float scale, out float bias)
        {
            float maxShadowDistance = shadowSettings.maxShadowDistance.value;
            float maxShadowDistanceSq = maxShadowDistance * maxShadowDistance;
            float cascadeBorder;
            int splitCount = shadowSettings.cascadeShadowSplitCount.value;
            if (splitCount == 4)
                cascadeBorder = shadowSettings.cascadeShadowBorder3.value;
            else if (splitCount == 3)
                cascadeBorder = shadowSettings.cascadeShadowBorder2.value;
            else if (splitCount == 2)
                cascadeBorder = shadowSettings.cascadeShadowBorder1.value;
            else
                cascadeBorder = shadowSettings.cascadeShadowBorder0.value;

            GetScaleAndBiasForLinearDistanceFade(maxShadowDistanceSq, cascadeBorder, out scale, out bias);
        }

        void GetScaleAndBiasForLinearDistanceFade(float fadeDistance, float border, out float scale, out float bias)
        {
            // To avoid division from zero
            // This values ensure that fade within cascade will be 0 and outside 1
            if (border < 0.0001f)
            {
                float multiplier = 1000f; // To avoid blending if difference is in fractions
                scale = multiplier;
                bias = -fadeDistance * multiplier;
                return;
            }

            border = 1 - border;
            border *= border;

            // Fade with distance calculation is just a linear fade from 90% of fade distance to fade distance. 90% arbitrarily chosen but should work well enough.
            float distanceFadeNear = border * fadeDistance;
            scale = 1.0f / (fadeDistance - distanceFadeNear);
            bias = -distanceFadeNear / (fadeDistance - distanceFadeNear);
        }

        public void PushGlobalParameters(CommandBuffer cmd)
        {
            // This code must be in sync with HDShadowContext.hlsl
            cmd.SetGlobalBuffer(HDShaderIDs._HDShadowDatas, m_ShadowDataBuffer);
            cmd.SetGlobalBuffer(HDShaderIDs._HDDirectionalShadowData, m_DirectionalShadowDataBuffer);
        }

        public int GetShadowRequestCount()
        {
            return m_ShadowRequestCount;
        }

        public void Clear()
        {
            if (m_MaxShadowRequests == 0)
                return;

            // Clear the shadows atlas infos and requests
            m_Atlas.Clear();
            m_CascadeAtlas.Clear();
            if (ShaderConfig.s_AreaLights == 1)
                m_AreaLightShadowAtlas.Clear();

            cachedShadowManager.ClearShadowRequests();

            m_ShadowResolutionRequestCounter = 0;

            m_ShadowRequestCount = 0;
            m_CascadeCount = 0;
        }

        // Warning: must be called after ProcessShadowRequests and RenderShadows to have valid informations
        public void DisplayShadowAtlas(RTHandle atlasTexture, CommandBuffer cmd, Material debugMaterial, float screenX, float screenY, float screenSizeX, float screenSizeY, float minValue, float maxValue, MaterialPropertyBlock mpb)
        {
            m_Atlas.DisplayAtlas(atlasTexture, cmd, debugMaterial, new Rect(0, 0, m_Atlas.width, m_Atlas.height), screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
        }

        // Warning: must be called after ProcessShadowRequests and RenderShadows to have valid informations
        public void DisplayShadowCascadeAtlas(RTHandle atlasTexture, CommandBuffer cmd, Material debugMaterial, float screenX, float screenY, float screenSizeX, float screenSizeY, float minValue, float maxValue, MaterialPropertyBlock mpb)
        {
            m_CascadeAtlas.DisplayAtlas(atlasTexture, cmd, debugMaterial, new Rect(0, 0, m_CascadeAtlas.width, m_CascadeAtlas.height), screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
        }

        // Warning: must be called after ProcessShadowRequests and RenderShadows to have valid informations
        public void DisplayAreaLightShadowAtlas(RTHandle atlasTexture, CommandBuffer cmd, Material debugMaterial, float screenX, float screenY, float screenSizeX, float screenSizeY, float minValue, float maxValue, MaterialPropertyBlock mpb)
        {
            if (ShaderConfig.s_AreaLights == 1)
                m_AreaLightShadowAtlas.DisplayAtlas(atlasTexture, cmd, debugMaterial, new Rect(0, 0, m_AreaLightShadowAtlas.width, m_AreaLightShadowAtlas.height), screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
        }

        public void DisplayCachedPunctualShadowAtlas(RTHandle atlasTexture, CommandBuffer cmd, Material debugMaterial, float screenX, float screenY, float screenSizeX, float screenSizeY, float minValue, float maxValue, MaterialPropertyBlock mpb)
        {
            cachedShadowManager.punctualShadowAtlas.DisplayAtlas(atlasTexture, cmd, debugMaterial, new Rect(0, 0, cachedShadowManager.punctualShadowAtlas.width, cachedShadowManager.punctualShadowAtlas.height), screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
        }

        public void DisplayCachedAreaShadowAtlas(RTHandle atlasTexture, CommandBuffer cmd, Material debugMaterial, float screenX, float screenY, float screenSizeX, float screenSizeY, float minValue, float maxValue, MaterialPropertyBlock mpb)
        {
            if (ShaderConfig.s_AreaLights == 1)
                cachedShadowManager.areaShadowAtlas.DisplayAtlas(atlasTexture, cmd, debugMaterial, new Rect(0, 0, cachedShadowManager.areaShadowAtlas.width, cachedShadowManager.areaShadowAtlas.height), screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
        }

        // Warning: must be called after ProcessShadowRequests and RenderShadows to have valid informations
        public void DisplayShadowMap(in ShadowResult atlasTextures, int shadowIndex, CommandBuffer cmd, Material debugMaterial, float screenX, float screenY, float screenSizeX, float screenSizeY, float minValue, float maxValue, MaterialPropertyBlock mpb)
        {
            if (shadowIndex >= m_ShadowRequestCount)
                return;

            HDShadowRequest shadowRequest = m_ShadowRequests[shadowIndex];

            switch (shadowRequest.shadowMapType)
            {
                case ShadowMapType.PunctualAtlas:
                {
                    if (shadowRequest.isInCachedAtlas)
                        cachedShadowManager.punctualShadowAtlas.DisplayAtlas(atlasTextures.cachedPunctualShadowResult, cmd, debugMaterial, shadowRequest.cachedAtlasViewport, screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
                    else
                        m_Atlas.DisplayAtlas(atlasTextures.punctualShadowResult, cmd, debugMaterial, shadowRequest.dynamicAtlasViewport, screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
                    break;
                }
                case ShadowMapType.CascadedDirectional:
                {
                    m_CascadeAtlas.DisplayAtlas(atlasTextures.directionalShadowResult, cmd, debugMaterial, shadowRequest.dynamicAtlasViewport, screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
                    break;
                }
                case ShadowMapType.AreaLightAtlas:
                {
                    if (ShaderConfig.s_AreaLights == 1)
                    {
                        if (shadowRequest.isInCachedAtlas)
                            cachedShadowManager.areaShadowAtlas.DisplayAtlas(atlasTextures.cachedAreaShadowResult, cmd, debugMaterial, shadowRequest.cachedAtlasViewport, screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
                        else
                            m_AreaLightShadowAtlas.DisplayAtlas(atlasTextures.areaShadowResult, cmd, debugMaterial, shadowRequest.dynamicAtlasViewport, screenX, screenY, screenSizeX, screenSizeY, minValue, maxValue, mpb);
                    }
                    break;
                }
            }
            ;
        }

        internal static ShadowResult ReadShadowResult(in ShadowResult shadowResult, RenderGraphBuilder builder)
        {
            var result = new ShadowResult();

            if (shadowResult.punctualShadowResult.IsValid())
                result.punctualShadowResult = builder.ReadTexture(shadowResult.punctualShadowResult);
            if (shadowResult.directionalShadowResult.IsValid())
                result.directionalShadowResult = builder.ReadTexture(shadowResult.directionalShadowResult);
            if (shadowResult.areaShadowResult.IsValid())
                result.areaShadowResult = builder.ReadTexture(shadowResult.areaShadowResult);
            if (shadowResult.cachedPunctualShadowResult.IsValid())
                result.cachedPunctualShadowResult = builder.ReadTexture(shadowResult.cachedPunctualShadowResult);
            if (shadowResult.cachedAreaShadowResult.IsValid())
                result.cachedAreaShadowResult = builder.ReadTexture(shadowResult.cachedAreaShadowResult);

            return result;
        }

        internal void RenderShadows(RenderGraph renderGraph, in ShaderVariablesGlobal globalCB, HDCamera hdCamera, CullingResults cullResults, ref ShadowResult result)
        {
            InvalidateAtlasOutputsIfNeeded();

            // Avoid to do any commands if there is no shadow to draw
            if (m_ShadowRequestCount != 0 &&
                (hdCamera.frameSettings.IsEnabled(FrameSettingsField.OpaqueObjects) || hdCamera.frameSettings.IsEnabled(FrameSettingsField.TransparentObjects)))
            {
                // Punctual
                result.cachedPunctualShadowResult = cachedShadowManager.punctualShadowAtlas.RenderShadows(renderGraph, cullResults, globalCB, hdCamera.frameSettings, "Cached Punctual Lights Shadows rendering");
                BlitCachedShadows(renderGraph, ShadowMapType.PunctualAtlas);
                result.punctualShadowResult = m_Atlas.RenderShadows(renderGraph, cullResults, globalCB, hdCamera.frameSettings, "Punctual Lights Shadows rendering");

                if (ShaderConfig.s_AreaLights == 1)
                {
                    cachedShadowManager.areaShadowAtlas.RenderShadowMaps(renderGraph, cullResults, globalCB, hdCamera.frameSettings, "Cached Area Lights Shadows rendering");
                    BlitCachedShadows(renderGraph, ShadowMapType.AreaLightAtlas);
                    m_AreaLightShadowAtlas.RenderShadowMaps(renderGraph, cullResults, globalCB, hdCamera.frameSettings, "Area Light Shadows rendering");
                    result.areaShadowResult = m_AreaLightShadowAtlas.BlurShadows(renderGraph);
                    result.cachedAreaShadowResult = cachedShadowManager.areaShadowAtlas.BlurShadows(renderGraph);
                }

                if (cachedShadowManager.DirectionalHasCachedAtlas())
                {

                    if (cachedShadowManager.directionalLightAtlas.HasShadowRequests())
                    {
                        cachedShadowManager.UpdateDirectionalCacheTexture(renderGraph);
                        cachedShadowManager.directionalLightAtlas.RenderShadows(renderGraph, cullResults, globalCB, hdCamera.frameSettings, "Cached Directional Lights Shadows rendering");
                    }
                    BlitCachedShadows(renderGraph, ShadowMapType.CascadedDirectional);
                }
                result.directionalShadowResult = m_CascadeAtlas.RenderShadows(renderGraph, cullResults, globalCB, hdCamera.frameSettings, "Directional Light Shadows rendering");
            }

            // TODO RENDERGRAPH
            // Not really good to bind things globally here (makes lifecycle of the textures fuzzy)
            // Probably better to bind it explicitly where needed (deferred lighting and forward/debug passes)
            // We can probably remove this when we have only one code path and can clean things up a bit.
            BindShadowGlobalResources(renderGraph, result);
        }

        internal void ReleaseSharedShadowAtlases(RenderGraph renderGraph)
        {
            if (cachedShadowManager.DirectionalHasCachedAtlas())
                cachedShadowManager.directionalLightAtlas.CleanupRenderGraphOutput(renderGraph);

            cachedShadowManager.punctualShadowAtlas.CleanupRenderGraphOutput(renderGraph);
            if (ShaderConfig.s_AreaLights == 1)
                cachedShadowManager.areaShadowAtlas.CleanupRenderGraphOutput(renderGraph);

            cachedShadowManager.DefragAtlas(HDLightType.Point);
            cachedShadowManager.DefragAtlas(HDLightType.Spot);
            if (ShaderConfig.s_AreaLights == 1)
                cachedShadowManager.DefragAtlas(HDLightType.Area);
        }

        void InvalidateAtlasOutputsIfNeeded()
        {
            cachedShadowManager.punctualShadowAtlas.InvalidateOutputIfNeeded();
            m_Atlas.InvalidateOutputIfNeeded();
            m_CascadeAtlas.InvalidateOutputIfNeeded();
            if (cachedShadowManager.DirectionalHasCachedAtlas())
            {
                cachedShadowManager.directionalLightAtlas.InvalidateOutputIfNeeded();
            }
            if (ShaderConfig.s_AreaLights == 1)
            {
                cachedShadowManager.areaShadowAtlas.InvalidateOutputIfNeeded();
                m_AreaLightShadowAtlas.InvalidateOutputIfNeeded();
            }
        }

        class BindShadowGlobalResourcesPassData
        {
            public ShadowResult shadowResult;
        }


        static void BindAtlasTexture(RenderGraphContext ctx, TextureHandle texture, int shaderId)
        {
            if (texture.IsValid())
                ctx.cmd.SetGlobalTexture(shaderId, texture);
            else
                ctx.cmd.SetGlobalTexture(shaderId, ctx.defaultResources.defaultShadowTexture);
        }

        void BindShadowGlobalResources(RenderGraph renderGraph, in ShadowResult shadowResult)
        {
            using (var builder = renderGraph.AddRenderPass<BindShadowGlobalResourcesPassData>("BindShadowGlobalResources", out var passData))
            {
                passData.shadowResult = ReadShadowResult(shadowResult, builder);
                builder.AllowPassCulling(false);
                builder.SetRenderFunc(
                    (BindShadowGlobalResourcesPassData data, RenderGraphContext ctx) =>
                    {
                        BindAtlasTexture(ctx, data.shadowResult.punctualShadowResult, HDShaderIDs._ShadowmapAtlas);
                        BindAtlasTexture(ctx, data.shadowResult.directionalShadowResult, HDShaderIDs._ShadowmapCascadeAtlas);
                        BindAtlasTexture(ctx, data.shadowResult.areaShadowResult, HDShaderIDs._ShadowmapAreaAtlas);
                        BindAtlasTexture(ctx, data.shadowResult.cachedPunctualShadowResult, HDShaderIDs._CachedShadowmapAtlas);
                        BindAtlasTexture(ctx, data.shadowResult.cachedAreaShadowResult, HDShaderIDs._CachedAreaLightShadowmapAtlas);
                    });
            }
        }

        internal static void BindDefaultShadowGlobalResources(RenderGraph renderGraph)
        {
            using (var builder = renderGraph.AddRenderPass<BindShadowGlobalResourcesPassData>("BindDefaultShadowGlobalResources", out var passData))
            {
                builder.AllowPassCulling(false);
                builder.SetRenderFunc(
                    (BindShadowGlobalResourcesPassData data, RenderGraphContext ctx) =>
                    {
                        BindAtlasTexture(ctx, ctx.defaultResources.defaultShadowTexture, HDShaderIDs._ShadowmapAtlas);
                        BindAtlasTexture(ctx, ctx.defaultResources.defaultShadowTexture, HDShaderIDs._ShadowmapCascadeAtlas);
                        BindAtlasTexture(ctx, ctx.defaultResources.defaultShadowTexture, HDShaderIDs._ShadowmapAreaAtlas);
                        BindAtlasTexture(ctx, ctx.defaultResources.defaultShadowTexture, HDShaderIDs._CachedShadowmapAtlas);
                        BindAtlasTexture(ctx, ctx.defaultResources.defaultShadowTexture, HDShaderIDs._CachedAreaLightShadowmapAtlas);
                    });
            }
        }

        void BlitCachedShadows(RenderGraph renderGraph)
        {
            m_Atlas.BlitCachedIntoAtlas(renderGraph, cachedShadowManager.punctualShadowAtlas.GetOutputTexture(renderGraph), new Vector2Int(cachedShadowManager.punctualShadowAtlas.width, cachedShadowManager.punctualShadowAtlas.height), m_BlitShadowMaterial, "Blit Punctual Mixed Cached Shadows", HDProfileId.BlitPunctualMixedCachedShadowMaps);
            if (cachedShadowManager.DirectionalHasCachedAtlas())
            {
                m_CascadeAtlas.BlitCachedIntoAtlas(renderGraph, cachedShadowManager.directionalLightAtlas.GetOutputTexture(renderGraph), new Vector2Int(cachedShadowManager.directionalLightAtlas.width, cachedShadowManager.directionalLightAtlas.height), m_BlitShadowMaterial, "Blit Directional Mixed Cached Shadows", HDProfileId.BlitDirectionalMixedCachedShadowMaps);
            }

            if (ShaderConfig.s_AreaLights == 1)
            {
                m_AreaLightShadowAtlas.BlitCachedIntoAtlas(renderGraph, cachedShadowManager.areaShadowAtlas.GetOutputTexture(renderGraph), new Vector2Int(cachedShadowManager.areaShadowAtlas.width, cachedShadowManager.areaShadowAtlas.height), m_BlitShadowMaterial, "Blit Area Mixed Cached Shadows", HDProfileId.BlitAreaMixedCachedShadowMaps);
            }
        }

        void BlitCachedShadows(RenderGraph renderGraph, ShadowMapType shadowAtlas)
        {
            if (shadowAtlas == ShadowMapType.PunctualAtlas)
                m_Atlas.BlitCachedIntoAtlas(renderGraph, cachedShadowManager.punctualShadowAtlas.GetOutputTexture(renderGraph), new Vector2Int(cachedShadowManager.punctualShadowAtlas.width, cachedShadowManager.punctualShadowAtlas.height), m_BlitShadowMaterial, "Blit Punctual Mixed Cached Shadows", HDProfileId.BlitPunctualMixedCachedShadowMaps);
            if (shadowAtlas == ShadowMapType.CascadedDirectional && cachedShadowManager.DirectionalHasCachedAtlas())
                m_CascadeAtlas.BlitCachedIntoAtlas(renderGraph, cachedShadowManager.directionalLightAtlas.GetOutputTexture(renderGraph), new Vector2Int(cachedShadowManager.directionalLightAtlas.width, cachedShadowManager.directionalLightAtlas.height), m_BlitShadowMaterial, "Blit Directional Mixed Cached Shadows", HDProfileId.BlitDirectionalMixedCachedShadowMaps);
            if (shadowAtlas == ShadowMapType.AreaLightAtlas && ShaderConfig.s_AreaLights == 1)
                m_AreaLightShadowAtlas.BlitCachedIntoAtlas(renderGraph, cachedShadowManager.areaShadowAtlas.GetShadowMapDepthTexture(renderGraph), new Vector2Int(cachedShadowManager.areaShadowAtlas.width, cachedShadowManager.areaShadowAtlas.height), m_BlitShadowMaterial, "Blit Area Mixed Cached Shadows", HDProfileId.BlitAreaMixedCachedShadowMaps);
        }
    }
}
