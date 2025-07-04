using System;
using System.Reflection;
using System.Collections.Generic;

namespace UnityEngine.Rendering.HighDefinition
{
    [HDRPHelpURL("Default-Settings-Window")]
    partial class HDRenderPipelineRuntimeResources : HDRenderPipelineResources
    {
        [Serializable, ReloadGroup]
        public sealed class ShaderResources
        {
            // Defaults
            [Reload("Runtime/Material/Lit/Lit.shader")]
            public Shader defaultPS;

            // Debug
            [Reload("Runtime/Debug/DebugDisplayLatlong.Shader")]
            public Shader debugDisplayLatlongPS;
            [Reload("Runtime/Debug/DebugViewMaterialGBuffer.Shader")]
            public Shader debugViewMaterialGBufferPS;
            [Reload("Runtime/Debug/DebugViewTiles.Shader")]
            public Shader debugViewTilesPS;
            [Reload("Runtime/Debug/DebugFullScreen.Shader")]
            public Shader debugFullScreenPS;
            [Reload("Runtime/Debug/DebugColorPicker.Shader")]
            public Shader debugColorPickerPS;
            [Reload("Runtime/Debug/DebugExposure.Shader")]
            public Shader debugExposurePS;
            [Reload("Runtime/Debug/DebugHDR.Shader")]
            public Shader debugHDRPS;
            [Reload("Runtime/Debug/DebugLightVolumes.Shader")]
            public Shader debugLightVolumePS;
            [Reload("Runtime/Debug/DebugLightVolumes.compute")]
            public ComputeShader debugLightVolumeCS;
            [Reload("Runtime/Debug/DebugBlitQuad.Shader")]
            public Shader debugBlitQuad;
            [Reload("Runtime/Debug/DebugVTBlit.Shader")]
            public Shader debugViewVirtualTexturingBlit;
            [Reload("Runtime/Debug/MaterialError.Shader")]
            public Shader materialError;
            [Reload("Runtime/Debug/MaterialLoading.shader")]
            public Shader materialLoading;
            [Reload("Runtime/Debug/ClearDebugBuffer.compute")]
            public ComputeShader clearDebugBufferCS;

            // APV
            [Reload("Runtime/Debug/ProbeVolumeDebug.shader")]
            public Shader probeVolumeDebugShader;
            [Reload("Runtime/Debug/ProbeVolumeOffsetDebug.shader")]
            public Shader probeVolumeOffsetDebugShader;
            [Reload("Runtime/Lighting/ProbeVolume/ProbeVolumeBlendStates.compute")]
            public ComputeShader probeVolumeBlendStatesCS;

            [Reload("Runtime/Debug/DebugWaveform.shader")]
            public Shader debugWaveformPS;
            [Reload("Runtime/Debug/DebugWaveform.compute")]
            public ComputeShader debugWaveformCS;

            [Reload("Runtime/Debug/DebugVectorscope.shader")]
            public Shader debugVectorscopePS;
            [Reload("Runtime/Debug/DebugVectorscope.compute")]
            public ComputeShader debugVectorscopeCS;

            // Lighting
            [Reload("Runtime/Lighting/Deferred.Shader")]
            public Shader deferredPS;
            [Reload("Runtime/RenderPipeline/RenderPass/ColorPyramidPS.Shader")]
            public Shader colorPyramidPS;
            [Reload("Runtime/RenderPipeline/RenderPass/DepthPyramid.compute")]
            public ComputeShader depthPyramidCS;
            [Reload("Runtime/RenderPipeline/RenderPass/GenerateMaxZ.compute")]
            public ComputeShader maxZCS;

            [Reload("Runtime/Core/CoreResources/GPUCopy.compute")]
            public ComputeShader copyChannelCS;
            [Reload("Runtime/Lighting/ScreenSpaceLighting/ScreenSpaceReflections.compute")]
            public ComputeShader screenSpaceReflectionsCS;
            [Reload("Runtime/RenderPipeline/RenderPass/Distortion/ApplyDistortion.shader")]
            public Shader applyDistortionPS;

            // Lighting tile pass
            [Reload("Runtime/Lighting/LightLoop/cleardispatchindirect.compute")]
            public ComputeShader clearDispatchIndirectCS;
            [Reload("Runtime/Lighting/LightLoop/ClearLightLists.compute")]
            public ComputeShader clearLightListsCS;
            [Reload("Runtime/Lighting/LightLoop/builddispatchindirect.compute")]
            public ComputeShader buildDispatchIndirectCS;
            [Reload("Runtime/Lighting/LightLoop/scrbound.compute")]
            public ComputeShader buildScreenAABBCS;
            [Reload("Runtime/Lighting/LightLoop/lightlistbuild.compute")]
            public ComputeShader buildPerTileLightListCS;               // FPTL
            [Reload("Runtime/Lighting/LightLoop/lightlistbuild-bigtile.compute")]
            public ComputeShader buildPerBigTileLightListCS;
            [Reload("Runtime/Lighting/LightLoop/lightlistbuild-clustered.compute")]
            public ComputeShader buildPerVoxelLightListCS;              // clustered
            [Reload("Runtime/Lighting/LightLoop/lightlistbuild-clearatomic.compute")]
            public ComputeShader lightListClusterClearAtomicIndexCS;
            [Reload("Runtime/Lighting/LightLoop/materialflags.compute")]
            public ComputeShader buildMaterialFlagsCS;
            [Reload("Runtime/Lighting/LightLoop/Deferred.compute")]
            public ComputeShader deferredCS;
            [Reload("Runtime/Lighting/Shadow/ContactShadows.compute")]
            public ComputeShader contactShadowCS;
            [Reload("Runtime/Lighting/VolumetricLighting/VolumeVoxelization.compute")]
            public ComputeShader volumeVoxelizationCS;
            [Reload("Runtime/Lighting/VolumetricLighting/VolumetricLighting.compute")]
            public ComputeShader volumetricLightingCS;
            [Reload("Runtime/Lighting/VolumetricLighting/VolumetricLightingFiltering.compute")]
            public ComputeShader volumetricLightingFilteringCS;

            [Reload("Runtime/Lighting/LightLoop/DeferredTile.shader")]
            public Shader deferredTilePS;
            [Reload("Runtime/Lighting/Shadow/ScreenSpaceShadows.shader")]
            public Shader screenSpaceShadowPS;

            [Reload("Runtime/Material/SubsurfaceScattering/SubsurfaceScattering.compute")]
            public ComputeShader subsurfaceScatteringCS;                // Disney SSS
            [Reload("Runtime/Material/SubsurfaceScattering/CombineLighting.shader")]
            public Shader combineLightingPS;

            [Reload("Runtime/Lighting/VolumetricLighting/DebugLocalVolumetricFogAtlas.shader")]
            public Shader debugLocalVolumetricFogAtlasPS;

            // General
            [Reload("Runtime/RenderPipeline/RenderPass/MotionVectors/CameraMotionVectors.shader")]
            public Shader cameraMotionVectorsPS;
            [Reload("Runtime/ShaderLibrary/ClearStencilBuffer.shader")]
            public Shader clearStencilBufferPS;
            [Reload("Runtime/ShaderLibrary/CopyStencilBuffer.shader")]
            public Shader copyStencilBufferPS;
            [Reload("Runtime/ShaderLibrary/CopyDepthBuffer.shader")]
            public Shader copyDepthBufferPS;
            [Reload("Runtime/ShaderLibrary/Blit.shader")]
            public Shader blitPS;
            [Reload("Runtime/ShaderLibrary/BlitColorAndDepth.shader")]
            public Shader blitColorAndDepthPS;

            [Reload("Runtime/Core/CoreResources/ClearBuffer2D.compute")]
            public ComputeShader clearBuffer2D;

            [Reload("Runtime/ShaderLibrary/DownsampleDepth.shader")]
            public Shader downsampleDepthPS;
            [Reload("Runtime/ShaderLibrary/UpsampleTransparent.shader")]
            public Shader upsampleTransparentPS;

            [Reload("Runtime/ShaderLibrary/ResolveStencilBuffer.compute")]
            public ComputeShader resolveStencilCS;

            // Sky
            [Reload("Runtime/Sky/BlitCubemap.shader")]
            public Shader blitCubemapPS;
            [Reload("Runtime/Material/GGXConvolution/BuildProbabilityTables.compute")]
            public ComputeShader buildProbabilityTablesCS;
            [Reload("Runtime/Material/GGXConvolution/ComputeGgxIblSampleData.compute")]
            public ComputeShader computeGgxIblSampleDataCS;
            [Reload("Runtime/Material/GGXConvolution/GGXConvolve.shader")]
            public Shader GGXConvolvePS;
            [Reload("Runtime/Material/Fabric/CharlieConvolve.shader")]
            public Shader charlieConvolvePS;
            [Reload("Runtime/Lighting/AtmosphericScattering/OpaqueAtmosphericScattering.shader")]
            public Shader opaqueAtmosphericScatteringPS;
            [Reload("Runtime/Material/Glints/RadianceLevels/ComputeRadianceLevelWeights.shader")]
            public Shader computeRadianceLevelWeightsPS;
            [Reload("Runtime/Material/Glints/RadianceLevels/ComputeRadianceLevels.compute")]
            public ComputeShader computeRadianceLevelsCS;
            [Reload("Runtime/Sky/HDRISky/HDRISky.shader")]
            public Shader hdriSkyPS;
            [Reload("Runtime/Sky/HDRISky/IntegrateHDRISky.shader")]
            public Shader integrateHdriSkyPS;
            [Reload("Skybox/Cubemap", ReloadAttribute.Package.Builtin)]
            public Shader skyboxCubemapPS;
            [Reload("Runtime/Sky/GradientSky/GradientSky.shader")]
            public Shader gradientSkyPS;
            [Reload("Runtime/Sky/AmbientProbeConvolution.compute")]
            public ComputeShader ambientProbeConvolutionCS;
            [Reload("Runtime/Sky/PhysicallyBasedSky/SkyLUTGenerator.compute")]
            public ComputeShader skyLUTGenerator;
            [Reload("Runtime/Sky/PhysicallyBasedSky/GroundIrradiancePrecomputation.compute")]
            public ComputeShader groundIrradiancePrecomputationCS;
            [Reload("Runtime/Sky/PhysicallyBasedSky/InScatteredRadiancePrecomputation.compute")]
            public ComputeShader inScatteredRadiancePrecomputationCS;
            [Reload("Runtime/Sky/PhysicallyBasedSky/PhysicallyBasedSky.shader")]
            public Shader physicallyBasedSkyPS;
            [Reload("Runtime/Lighting/PlanarReflectionFiltering.compute")]
            public ComputeShader planarReflectionFilteringCS;
            [Reload("Runtime/Sky/CloudSystem/CloudLayer/CloudLayer.shader")]
            public Shader cloudLayerPS;
            [Reload("Runtime/Sky/CloudSystem/CloudLayer/BakeCloudTexture.compute")]
            public ComputeShader bakeCloudTextureCS;
            [Reload("Runtime/Sky/CloudSystem/CloudLayer/BakeCloudShadows.compute")]
            public ComputeShader bakeCloudShadowsCS;

            // Volumetric Clouds
            [Reload("Runtime/Lighting/VolumetricLighting/VolumetricClouds.compute")]
            public ComputeShader volumetricCloudsCS;
            [Reload("Editor/Lighting/VolumetricClouds/CloudMapGenerator.compute")]
            public ComputeShader volumetricCloudMapGeneratorCS;
            [Reload("Runtime/Lighting/VolumetricLighting/VolumetricCloudsCombine.shader")]
            public Shader volumetricCloudsCombinePS;

            // Water
            [Reload("Runtime/Water/WaterSimulation.compute")]
            public ComputeShader waterSimulationCS;
            [Reload("Runtime/Water/FourierTransform.compute")]
            public ComputeShader fourierTransformCS;
            [Reload("Runtime/RenderPipelineResources/ShaderGraph/Water.shadergraph")]
            public Shader waterPS;
            [Reload("Runtime/Water/WaterLighting.compute")]
            public ComputeShader waterLightingCS;
            [Reload("Runtime/Water/WaterCaustics.shader")]
            public Shader waterCausticsPS;

            // Material
            [Reload("Runtime/Material/PreIntegratedFGD/PreIntegratedFGD_GGXDisneyDiffuse.shader")]
            public Shader preIntegratedFGD_GGXDisneyDiffusePS;
            [Reload("Runtime/Material/PreIntegratedFGD/PreIntegratedFGD_CharlieFabricLambert.shader")]
            public Shader preIntegratedFGD_CharlieFabricLambertPS;
            [Reload("Runtime/Material/AxF/PreIntegratedFGD_Ward.shader")]
            public Shader preIntegratedFGD_WardPS;
            [Reload("Runtime/Material/AxF/PreIntegratedFGD_CookTorrance.shader")]
            public Shader preIntegratedFGD_CookTorrancePS;
            [Reload("Runtime/Material/PreIntegratedFGD/PreIntegratedFGD_Marschner.shader")]
            public Shader preIntegratedFGD_MarschnerPS;
            [Reload("Runtime/Material/PreIntegratedFGD/PreIntegratedFGD_DGGXOnly.shader")]
            public Shader preIntegratedFGD_DGGXOnlyPS;
            [Reload("Runtime/Material/Hair/MultipleScattering/HairMultipleScatteringPreIntegration.compute")]
            public ComputeShader preIntegratedFiberScatteringCS;
            [Reload("Runtime/Material/VolumetricMaterial/VolumetricMaterial.compute")]
            public ComputeShader volumetricMaterialCS;
            // Utilities / Core
            [Reload("Runtime/Core/CoreResources/EncodeBC6H.compute")]
            public ComputeShader encodeBC6HCS;
            [Reload("Runtime/Core/CoreResources/CubeToPano.shader")]
            public Shader cubeToPanoPS;
            [Reload("Runtime/Core/CoreResources/BlitCubeTextureFace.shader")]
            public Shader blitCubeTextureFacePS;
            [Reload("Runtime/Material/LTCAreaLight/FilterAreaLightCookies.shader")]
            public Shader filterAreaLightCookiesPS;
            [Reload("Runtime/Core/CoreResources/ClearUIntTextureArray.compute")]
            public ComputeShader clearUIntTextureCS;
            [Reload("Runtime/RenderPipeline/RenderPass/CustomPass/CustomPassUtils.shader")]
            public Shader customPassUtils;
            [Reload("Runtime/RenderPipeline/RenderPass/CustomPass/CustomPassRenderersUtils.shader")]
            public Shader customPassRenderersUtils;
            [Reload("Runtime/RenderPipeline/Utility/Texture3DAtlas.compute")]
            public ComputeShader texture3DAtlasCS;

            // XR
            [Reload("Runtime/ShaderLibrary/XRMirrorView.shader")]
            public Shader xrMirrorViewPS;
            [Reload("Runtime/ShaderLibrary/XROcclusionMesh.shader")]
            public Shader xrOcclusionMeshPS;

            // Shadow
            [Reload("Runtime/Lighting/Shadow/ShadowClear.shader")]
            public Shader shadowClearPS;
            [Reload("Runtime/Lighting/Shadow/EVSMBlur.compute")]
            public ComputeShader evsmBlurCS;
            [Reload("Runtime/Lighting/Shadow/DebugDisplayHDShadowMap.shader")]
            public Shader debugHDShadowMapPS;
            [Reload("Runtime/Lighting/Shadow/MomentShadows.compute")]
            public ComputeShader momentShadowsCS;
            [Reload("Runtime/Lighting/Shadow/ShadowBlit.shader")]
            public Shader shadowBlitPS;

            // Decal
            [Reload("Runtime/Material/Decal/DecalNormalBuffer.shader")]
            public Shader decalNormalBufferPS;

            // Ambient occlusion
            [Reload("Runtime/Lighting/ScreenSpaceLighting/GTAO.compute")]
            public ComputeShader GTAOCS;
            [Reload("Runtime/Lighting/ScreenSpaceLighting/GTAOSpatialDenoise.compute")]
            public ComputeShader GTAOSpatialDenoiseCS;
            [Reload("Runtime/Lighting/ScreenSpaceLighting/GTAOTemporalDenoise.compute")]
            public ComputeShader GTAOTemporalDenoiseCS;
            [Reload("Runtime/Lighting/ScreenSpaceLighting/GTAOCopyHistory.compute")]
            public ComputeShader GTAOCopyHistoryCS;
            [Reload("Runtime/Lighting/ScreenSpaceLighting/GTAOBlurAndUpsample.compute")]
            public ComputeShader GTAOBlurAndUpsample;

            [Reload("Runtime/Lighting/ScreenSpaceLighting/ScreenSpaceGlobalIllumination.compute")]
            public ComputeShader screenSpaceGlobalIlluminationCS;

            // MSAA Shaders
            [Reload("Runtime/RenderPipeline/RenderPass/MSAA/DepthValues.shader")]
            public Shader depthValuesPS;
            [Reload("Runtime/RenderPipeline/RenderPass/MSAA/ColorResolve.shader")]
            public Shader colorResolvePS;
            [Reload("Runtime/RenderPipeline/RenderPass/MSAA/MotionVecResolve.shader")]
            public Shader resolveMotionVecPS;

            // Post-processing
            [Reload("Runtime/PostProcessing/Shaders/AlphaCopy.compute")]
            public ComputeShader copyAlphaCS;
            [Reload("Runtime/PostProcessing/Shaders/NaNKiller.compute")]
            public ComputeShader nanKillerCS;
            [Reload("Runtime/PostProcessing/Shaders/Exposure.compute")]
            public ComputeShader exposureCS;
            [Reload("Runtime/PostProcessing/Shaders/HistogramExposure.compute")]
            public ComputeShader histogramExposureCS;
            [Reload("Runtime/PostProcessing/Shaders/ApplyExposure.compute")]
            public ComputeShader applyExposureCS;
            [Reload("Runtime/PostProcessing/Shaders/DebugHistogramImage.compute")]
            public ComputeShader debugImageHistogramCS;
            [Reload("Runtime/PostProcessing/Shaders/DebugHDRxyMapping.compute")]
            public ComputeShader debugHDRxyMappingCS;
            [Reload("Runtime/PostProcessing/Shaders/UberPost.compute")]
            public ComputeShader uberPostCS;
            [Reload("Runtime/PostProcessing/Shaders/LutBuilder3D.compute")]
            public ComputeShader lutBuilder3DCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldKernel.compute")]
            public ComputeShader depthOfFieldKernelCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldCoC.compute")]
            public ComputeShader depthOfFieldCoCCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldCoCReproject.compute")]
            public ComputeShader depthOfFieldCoCReprojectCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldCoCDilate.compute")]
            public ComputeShader depthOfFieldDilateCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldMip.compute")]
            public ComputeShader depthOfFieldMipCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldMipSafe.compute")]
            public ComputeShader depthOfFieldMipSafeCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldPrefilter.compute")]
            public ComputeShader depthOfFieldPrefilterCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldTileMax.compute")]
            public ComputeShader depthOfFieldTileMaxCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldGather.compute")]
            public ComputeShader depthOfFieldGatherCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldCombine.compute")]
            public ComputeShader depthOfFieldCombineCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldPreCombineFar.compute")]
            public ComputeShader depthOfFieldPreCombineFarCS;
            [Reload("Runtime/PostProcessing/Shaders/DepthOfFieldClearIndirectArgs.compute")]
            public ComputeShader depthOfFieldClearIndirectArgsCS;
            [Reload("Runtime/PostProcessing/Shaders/PaniniProjection.compute")]
            public ComputeShader paniniProjectionCS;
            [Reload("Runtime/PostProcessing/Shaders/MotionBlurMotionVecPrep.compute")]
            public ComputeShader motionBlurMotionVecPrepCS;
            [Reload("Runtime/PostProcessing/Shaders/MotionBlurGenTilePass.compute")]
            public ComputeShader motionBlurGenTileCS;
            [Reload("Runtime/PostProcessing/Shaders/MotionBlurMergeTilePass.compute")]
            public ComputeShader motionBlurMergeTileCS;
            [Reload("Runtime/PostProcessing/Shaders/MotionBlurNeighborhoodTilePass.compute")]
            public ComputeShader motionBlurNeighborhoodTileCS;
            [Reload("Runtime/PostProcessing/Shaders/MotionBlur.compute")]
            public ComputeShader motionBlurCS;
            [Reload("Runtime/PostProcessing/Shaders/BloomPrefilter.compute")]
            public ComputeShader bloomPrefilterCS;
            [Reload("Runtime/PostProcessing/Shaders/BloomBlur.compute")]
            public ComputeShader bloomBlurCS;
            [Reload("Runtime/PostProcessing/Shaders/BloomUpsample.compute")]
            public ComputeShader bloomUpsampleCS;
            [Reload("Runtime/PostProcessing/Shaders/FXAA.compute")]
            public ComputeShader FXAACS;
            [Reload("Runtime/PostProcessing/Shaders/FinalPass.shader")]
            public Shader finalPassPS;
            [Reload("Runtime/PostProcessing/Shaders/ClearBlack.shader")]
            public Shader clearBlackPS;
            [Reload("Runtime/PostProcessing/Shaders/SubpixelMorphologicalAntialiasing.shader")]
            public Shader SMAAPS;
            [Reload("Runtime/PostProcessing/Shaders/TemporalAntialiasing.shader")]
            public Shader temporalAntialiasingPS;
            [Reload("Runtime/PostProcessing/Shaders/LensFlareDataDriven.shader")]
            public Shader lensFlareDataDrivenPS;
            [Reload("Runtime/PostProcessing/Shaders/LensFlareMergeOcclusionDataDriven.compute")]
            public ComputeShader lensFlareMergeOcclusionCS;
            [Reload("Runtime/PostProcessing/Shaders/DLSSBiasColorMask.shader")]
            public Shader DLSSBiasColorMaskPS;
            [Reload("Runtime/PostProcessing/Shaders/CompositeWithUIAndOETF.shader")]
            public Shader compositeUIAndOETFApplyPS;

            // Physically based DoF
            [Reload("Runtime/PostProcessing/Shaders/DoFCircleOfConfusion.compute")]
            public ComputeShader dofCircleOfConfusion;
            [Reload("Runtime/PostProcessing/Shaders/DoFGather.compute")]
            public ComputeShader dofGatherCS;
            [Reload("Runtime/PostProcessing/Shaders/DoFCoCMinMax.compute")]
            public ComputeShader dofCoCMinMaxCS;
            [Reload("Runtime/PostProcessing/Shaders/DoFMinMaxDilate.compute")]
            public ComputeShader dofMinMaxDilateCS;
            [Reload("Runtime/PostProcessing/Shaders/DoFCombine.compute")]
            public ComputeShader dofCombineCS;

            [Reload("Runtime/PostProcessing/Shaders/ContrastAdaptiveSharpen.compute")]
            public ComputeShader contrastAdaptiveSharpenCS;
            [Reload("Runtime/PostProcessing/Shaders/EdgeAdaptiveSpatialUpsampling.compute")]
            public ComputeShader edgeAdaptiveSpatialUpsamplingCS;
            [Reload("Runtime/VirtualTexturing/Shaders/DownsampleVTFeedback.compute")]
            public ComputeShader VTFeedbackDownsample;

            // Accumulation
            [Reload("Runtime/RenderPipeline/Accumulation/Shaders/Accumulation.compute")]
            public ComputeShader accumulationCS;

            [Reload("Runtime/RenderPipeline/Accumulation/Shaders/BlitAndExpose.compute")]
            public ComputeShader blitAndExposeCS;

            // Compositor
            [Reload("Runtime/Compositor/Shaders/AlphaInjection.shader")]
            public Shader alphaInjectionPS;
            [Reload("Runtime/Compositor/Shaders/ChromaKeying.shader")]
            public Shader chromaKeyingPS;
            [Reload("Runtime/Compositor/Shaders/CustomClear.shader")]
            public Shader customClearPS;

            // Denoising
            [Reload("Runtime/Lighting/ScreenSpaceLighting/BilateralUpsample.compute")]
            public ComputeShader bilateralUpsampleCS;
            [Reload("Runtime/RenderPipeline/Raytracing/Shaders/Denoising/TemporalFilter.compute")]
            public ComputeShader temporalFilterCS;
            [Reload("Runtime/RenderPipeline/Raytracing/Shaders/Denoising/DiffuseDenoiser.compute")]
            public ComputeShader diffuseDenoiserCS;

#if UNITY_EDITOR
            // Furnace Testing (BSDF Energy Conservation)
            [Reload("Tests/Editor/Utilities/FurnaceTests.compute")]
            public ComputeShader furnaceTestCS;
#endif

#if UNITY_EDITOR
            // Iterator to retrieve all compute shaders in reflection so we don't have to keep a list of
            // used compute shaders up to date (prefer editor-only usage)
            public IEnumerable<ComputeShader> GetAllComputeShaders()
            {
                var fields = typeof(ShaderResources).GetFields(BindingFlags.Public | BindingFlags.Instance);

                foreach (var field in fields)
                {
                    if (field.GetValue(this) is ComputeShader computeShader)
                        yield return computeShader;
                }
            }

#endif
        }

        [Serializable, ReloadGroup]
        public sealed class MaterialResources
        {
			[Reload("Runtime/RenderPipelineResources/Material/AreaLightCookieViewer.mat")]
            public Material areaLightCookieMaterial; // We also need one for the cookie because the emissive map is a keyword in our Unlit shader.
        }

        [Serializable, ReloadGroup]
        public sealed class TextureResources
        {
            // Debug
            [Reload("Runtime/RenderPipelineResources/Texture/DebugFont.tga")]
            public Texture2D debugFontTex;
            [Reload("Runtime/Debug/ColorGradient.png")]
            public Texture2D colorGradient;
            [Reload("Runtime/RenderPipelineResources/Texture/Matcap/DefaultMatcap.png")]
            public Texture2D matcapTex;

            // Pre-baked noise
            [Reload("Runtime/RenderPipelineResources/Texture/BlueNoise16/L/LDR_LLL1_{0}.png", 0, 32)]
            public Texture2D[] blueNoise16LTex;
            [Reload("Runtime/RenderPipelineResources/Texture/BlueNoise16/RGB/LDR_RGB1_{0}.png", 0, 32)]
            public Texture2D[] blueNoise16RGBTex;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/OwenScrambledNoise4.png")]
            public Texture2D owenScrambledRGBATex;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/OwenScrambledNoise256.png")]
            public Texture2D owenScrambled256Tex;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/ScrambleNoise.png")]
            public Texture2D scramblingTex;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/RankingTile1SPP.png")]
            public Texture2D rankingTile1SPP;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/ScramblingTile1SPP.png")]
            public Texture2D scramblingTile1SPP;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/RankingTile8SPP.png")]
            public Texture2D rankingTile8SPP;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/ScramblingTile8SPP.png")]
            public Texture2D scramblingTile8SPP;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/RankingTile256SPP.png")]
            public Texture2D rankingTile256SPP;
            [Reload("Runtime/RenderPipelineResources/Texture/CoherentNoise/ScramblingTile256SPP.png")]
            public Texture2D scramblingTile256SPP;

            // Precalculated eye caustic LUT
            [Reload("Runtime/RenderPipelineResources/Texture/EyeCausticLUT16R.exr")]
            public Texture3D eyeCausticLUT;

            // Clouds textures
            [Reload("Runtime/RenderPipelineResources/Texture/VolumetricClouds/CloudLutRainAO.png")]
            public Texture2D cloudLutRainAO;
            [Reload("Runtime/RenderPipelineResources/Texture/VolumetricClouds/WorleyNoise128RGBA.png")]
            public Texture3D worleyNoise128RGBA;
            [Reload("Runtime/RenderPipelineResources/Texture/VolumetricClouds/WorleyNoise32RGB.png")]
            public Texture3D worleyNoise32RGB;
            [Reload("Runtime/RenderPipelineResources/Texture/VolumetricClouds/PerlinNoise32RGB.png")]
            public Texture3D perlinNoise32RGB;

            // Water textures
            [Reload("Runtime/RenderPipelineResources/Texture/Water/FoamSurface.png")]
            public Texture2D foamSurface;

            // Post-processing
            [Reload(new[]
            {
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Thin01.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Thin02.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Medium01.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Medium02.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Medium03.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Medium04.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Medium05.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Medium06.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Large01.png",
                "Runtime/RenderPipelineResources/Texture/FilmGrain/Large02.png"
            })]
            public Texture2D[] filmGrainTex;
            [Reload("Runtime/RenderPipelineResources/Texture/SMAA/SearchTex.tga")]
            public Texture2D SMAASearchTex;
            [Reload("Runtime/RenderPipelineResources/Texture/SMAA/AreaTex.tga")]
            public Texture2D SMAAAreaTex;

            [Reload("Runtime/RenderPipelineResources/Texture/DefaultHDRISky.exr")]
            public Cubemap defaultHDRISky;

            [Reload("Runtime/RenderPipelineResources/Texture/DefaultCloudMap.png")]
            public Texture2D defaultCloudMap;
        }

        [Serializable, ReloadGroup]
        public sealed class ShaderGraphResources
        {
            [Reload("Runtime/ShaderLibrary/SolidColor.shadergraph")]
            public Shader objectIDPS;
            [Reload("Runtime/RenderPipelineResources/ShaderGraph/DefaultFogVolume.shadergraph")]
            public Shader defaultFogVolumeShader;
        }

        [Serializable, ReloadGroup]
        public sealed class AssetResources
        {
            [Reload("Runtime/RenderPipelineResources/defaultDiffusionProfile.asset")]
            public DiffusionProfileSettings defaultDiffusionProfile;

            //Area Light Emissive Meshes
            [Reload("Runtime/RenderPipelineResources/Mesh/Cylinder.fbx")]
            public Mesh emissiveCylinderMesh;
            [Reload("Runtime/RenderPipelineResources/Mesh/Quad.fbx")]
            public Mesh emissiveQuadMesh;
            [Reload("Runtime/RenderPipelineResources/Mesh/Sphere.fbx")]
            public Mesh sphereMesh;
            [Reload("Runtime/RenderPipelineResources/Mesh/ProbeDebugSphere.fbx")]
            public Mesh probeDebugSphere;
            [Reload("Runtime/RenderPipelineResources/Mesh/ProbeDebugPyramid.fbx")]
            public Mesh pyramidMesh;
        }

        public ShaderResources shaders;
        public MaterialResources materials;
        public TextureResources textures;
        public ShaderGraphResources shaderGraphs;
        public AssetResources assets;
    }
}
