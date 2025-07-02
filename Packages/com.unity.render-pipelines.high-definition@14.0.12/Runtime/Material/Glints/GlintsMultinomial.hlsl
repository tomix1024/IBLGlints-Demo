
// Include Deliot & Belcour implementation
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Glints2023.hlsl"

void SampleBinomialValueForSurfaceCellMultinomial(float randB, float randG, float successProb, float microfacetCount, out float resultPos, out float resultNeg)
{
#if defined(_GLINTS_IBL_DISTRIBUTION_UNGATED)
	resultPos = EvaluateBinomialValueUngated(randG, successProb, microfacetCount);
	resultNeg = microfacetCount - resultPos;
#elif defined(_GLINTS_IBL_DISTRIBUTION_GATED)
	resultPos = EvaluateBinomialValueGated(randB, randG, successProb, microfacetCount);
	resultNeg = microfacetCount - resultPos;
#else
	EvaluateBinomialValueDualGated(randB, randG, successProb, microfacetCount, resultPos, resultNeg);
#endif
}

void SampleBinomialValueForSurfaceCellMultinomial(float2 randB, float2 randG, float2 successProb, float2 microfacetCount, out float2 resultPos, out float2 resultNeg)
{
#if defined(_GLINTS_IBL_DISTRIBUTION_UNGATED)
	resultPos = EvaluateBinomialValueUngated(randG, successProb, microfacetCount);
	resultNeg = microfacetCount - resultPos;
#elif defined(_GLINTS_IBL_DISTRIBUTION_GATED)
	resultPos = EvaluateBinomialValueGated(randB, randG, successProb, microfacetCount);
	resultNeg = microfacetCount - resultPos;
#else
	EvaluateBinomialValueDualGated(randB, randG, successProb, microfacetCount, resultPos, resultNeg);
#endif
}

void SampleBinomialValueForSurfaceCellMultinomial(float4 randB, float4 randG, float4 successProb, float4 microfacetCount, out float4 resultPos, out float4 resultNeg)
{
#if defined(_GLINTS_IBL_DISTRIBUTION_UNGATED)
	resultPos = EvaluateBinomialValueUngated(randG, successProb, microfacetCount);
	resultNeg = microfacetCount - resultPos;
#elif defined(_GLINTS_IBL_DISTRIBUTION_GATED)
	resultPos = EvaluateBinomialValueGated(randB, randG, successProb, microfacetCount);
	resultNeg = microfacetCount - resultPos;
#else
	EvaluateBinomialValueDualGated(randB, randG, successProb, microfacetCount, resultPos, resultNeg);
#endif
}


float4 SampleGlintGridSimplexVertexMultinomial(float4 randSlopesB, float4 randSlopesG, float4 successProbs, float microfacetCountBlended)
{
	// NOTE: Require successProbs to be normalized!
	float successProbA = successProbs.x + successProbs.y;
	float successProbB = max(0, 1 - successProbA);
	float successProb0inA = successProbs.x / max(1e-8, successProbA);
	float successProb2inB = successProbs.z / max(1e-8, successProbB);

	float result_A, result_B;
	SampleBinomialValueForSurfaceCellMultinomial(randSlopesB.x, randSlopesG.x, successProbA, microfacetCountBlended, result_A, result_B);
	//result_A = successProbA * microfacetCountBlended;
	//result_B = (1-successProbA) * microfacetCountBlended;

	// TODO these can be parallelized into float2's
	float result_0, result_1;
	SampleBinomialValueForSurfaceCellMultinomial(randSlopesB.y, randSlopesG.y, successProb0inA, result_A, result_0, result_1);
	float result_2, result_3;
	SampleBinomialValueForSurfaceCellMultinomial(randSlopesB.z, randSlopesG.z, successProb2inB, result_B, result_2, result_3);
	//result_0 = successProb0inA * result_A;
	//result_1 = (1-successProb0inA) * result_A;
	//result_2 = successProb2inB * result_B;
	//result_3 = (1-successProb2inB) * result_B;

	float4 result;
	result.x = result_0;
	result.y = result_1;
	result.z = result_2;
	result.w = result_3;
	return result;
}

float4 SampleGlintGridSimplexMultinomial(float2 uv, uint gridSeed, float footprintArea, float4 successProbs, float gridWeight)
{
	int2 glint0, glint1, glint2;
	float3 barycentrics;
	EvalSkewedGrid(uv, glint0, glint1, glint2, barycentrics);

	// Get per surface cell per slope cell random numbers
	uint3 seed0 = uint3(glint0 + 2147483648, gridSeed);
	uint3 seed1 = uint3(glint1 + 2147483648, gridSeed);
	uint3 seed2 = uint3(glint2 + 2147483648, gridSeed);

	// Random numbers for angular counting
	seed0 = pcg3d(seed0);
	seed1 = pcg3d(seed1);
	seed2 = pcg3d(seed2);
	float3 rand0 = UintToUniformFloat(seed0); // pcg3dFloat(uint3(glint0 + 2147483648, gridSeed)); // TODO : optimize away manual seeds
	float3 rand1 = UintToUniformFloat(seed1); // pcg3dFloat(uint3(glint1 + 2147483648, gridSeed));
	float3 rand2 = UintToUniformFloat(seed2); // pcg3dFloat(uint3(glint2 + 2147483648, gridSeed));

	// Random numbers for footprintMicrofacetCount
	float3 rand0x = pcg3dFloat(seed0);
	float3 rand1x = pcg3dFloat(seed1);
	float3 rand2x = pcg3dFloat(seed2);

	// Compute microfacet count with randomization
	float3 logDensityRand = clamp(SampleNormalDistribution(float3(rand0.x, rand1.x, rand2.x), _LogMicrofacetDensity, _DensityRandomization), 0.0, 50.0); // TODO : optimize sampleNormalDist
	float3 microfacetCountTotal = exp(logDensityRand);
	float3 microfacetCount = microfacetCountTotal * gridWeight;
	float3 microfacetCountBlended = ComputeMicrofacetCountInFootprint(rand0x, rand1x, rand2x, footprintArea, microfacetCount);
	float3 microfacetCountAverage = microfacetCount*footprintArea;

	// Generate per surface cell random numbers
	float4 rand0SlopesB, rand1SlopesB, rand2SlopesB;
	float4 rand0SlopesG, rand1SlopesG, rand2SlopesG;
	float2 slopeLerpDummy;
	CustomRand4Texture(float2(0,0), rand0.yz, rand0SlopesB, rand0SlopesG, slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand1.yz, rand1SlopesB, rand1SlopesG, slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand2.yz, rand2SlopesB, rand2SlopesG, slopeLerpDummy);

	// NOTE: Require successProbs to be normalized!
	float successProbA = successProbs.x + successProbs.y;
	float successProbB = max(0, 1 - successProbA);
	float successProb0inA = successProbs.x / max(1e-8, successProbA);
	float successProb2inB = successProbs.z / max(1e-8, successProbB);

	// Generate numbers of reflecting microfacets
	float4 result0 = SampleGlintGridSimplexVertexMultinomial(rand0SlopesB, rand0SlopesG, successProbs, microfacetCountBlended.x);
	float4 result1 = SampleGlintGridSimplexVertexMultinomial(rand1SlopesB, rand1SlopesG, successProbs, microfacetCountBlended.y);
	float4 result2 = SampleGlintGridSimplexVertexMultinomial(rand2SlopesB, rand2SlopesG, successProbs, microfacetCountBlended.z);

	// Interpolate result for glint grid cell
	float3 results0 = float3(result0.x, result1.x, result2.x);// / max(1e-4, microfacetCount.xyz);
	float3 results1 = float3(result0.y, result1.y, result2.y);// / max(1e-4, microfacetCount.xyz);
	float3 results2 = float3(result0.z, result1.z, result2.z);// / max(1e-4, microfacetCount.xyz);
	float3 results3 = float3(result0.w, result1.w, result2.w);// / max(1e-4, microfacetCount.xyz);

	// TODO NO GRID WEIGHT IN THIS DIVISION!
	results0 = results0 / (microfacetCountTotal*footprintArea);
	results1 = results1 / (microfacetCountTotal*footprintArea);
	results2 = results2 / (microfacetCountTotal*footprintArea);
	results3 = results3 / (microfacetCountTotal*footprintArea);

	float4 result;
	result.x = dot(results0, barycentrics);
	result.y = dot(results1, barycentrics);
	result.z = dot(results2, barycentrics);
	result.w = dot(results3, barycentrics);
	return result;
}

float4 SampleGlintsMultinomialNDF_Internal(float4 successProbs, float2 uv, float2 duvdx, float2 duvdy)
{
	// ACCURATE PIXEL FOOTPRINT ELLIPSE
	float2 ellipseMajorAxis, ellipseMinorAxis;
	float ellipseMajorLength, ellipseMinorLength;
	GetGradientEllipseNew(duvdx, duvdy, ellipseMajorAxis, ellipseMinorAxis, ellipseMajorLength, ellipseMinorLength);
	float ellipseRatio = ellipseMajorLength / ellipseMinorLength;

	// MANUAL LOD COMPENSATION
	float lod = log2(ellipseMinorLength * _ScreenSpaceScale);
	float lod0 = floor(lod); //lod >= 0.0 ? (int)(lod) : (int)(lod - 1.0);
	float lod1 = lod0 + 1;
	float divLod0 = exp2(lod0);
	float divLod1 = exp2(lod1);
	float lodLerp = frac(lod);
	float footprintAreaLOD0 = divLod0*divLod0;// = pow(exp2(lod0), 2.0);
	float footprintAreaLOD1 = divLod1*divLod1;// = pow(exp2(lod1), 2.0);

	// MANUAL ANISOTROPY RATIO COMPENSATION
	float ratio0 = max(pow(2.0, (int)log2(ellipseRatio)), 1.0);
	float ratio1 = ratio0 * 2.0;
	float ratioLerp = clamp(Remap(ellipseRatio, ratio0, ratio1, 0.0, 1.0), 0.0, 1.0);

	// MANUAL ANISOTROPY ROTATION COMPENSATION
	float2 v1 = float2(0.0, 1.0);
	float2 v2 = ellipseMajorAxis;
	float theta = atan2(v1.x * v2.y - v1.y * v2.x, v1.x * v2.x + v1.y * v2.y) * RAD2DEG;
	// Sanitize invalid theta
	theta = max(theta, 0);
	float thetaGrid = 90.0 / max(ratio0, 2.0);
	float thetaBin = (int)(theta / thetaGrid) * thetaGrid;
	thetaBin = thetaBin + (thetaGrid / 2.0);
	float thetaBin0 = theta < thetaBin ? thetaBin - thetaGrid / 2.0 : thetaBin;
	float thetaBinH = thetaBin0 + thetaGrid / 4.0;
	float thetaBin1 = thetaBin0 + thetaGrid / 2.0;
	float thetaBinLerp = Remap(theta, thetaBin0, thetaBin1, 0.0, 1.0);
	thetaBin0 = thetaBin0 <= 0.0 ? 180.0 + thetaBin0 : thetaBin0;

	// TETRAHEDRONIZATION OF ROTATION + RATIO + LOD GRID
	bool centerSpecialCase = (ratio0 == 1.0);
	float2 divLods = float2(divLod0, divLod1);
	float2 footprintAreas = float2(footprintAreaLOD0, footprintAreaLOD1);
	float2 ratios = float2(ratio0, ratio1);
	float4 thetaBins = float4(thetaBin0, thetaBinH, thetaBin1, 0.0); // added 0.0 for center singularity case
	float3 tetraA, tetraB, tetraC, tetraD;
	GetAnisoCorrectingGridTetrahedron(centerSpecialCase, thetaBinLerp, ratioLerp, lodLerp, tetraA, tetraB, tetraC, tetraD);
	if (centerSpecialCase == true) // Account for center singularity in barycentric computation
		thetaBinLerp = Remap01To(thetaBinLerp, 0.0, ratioLerp);
	float4 tetraBarycentricWeights = GetBarycentricWeightsTetrahedron(float3(thetaBinLerp, ratioLerp, lodLerp), tetraA, tetraB, tetraC, tetraD); // Compute barycentric coordinates within chosen tetrahedron

	// PREPARE NEEDED ROTATIONS
	tetraA.x *= 2; tetraB.x *= 2; tetraC.x *= 2; tetraD.x *= 2;
	if (centerSpecialCase == true) // Account for center singularity (if center vertex => no rotation)
	{
		tetraA.x = (tetraA.y == 0) ? 3 : tetraA.x;
		tetraB.x = (tetraB.y == 0) ? 3 : tetraB.x;
		tetraC.x = (tetraC.y == 0) ? 3 : tetraC.x;
		tetraD.x = (tetraD.y == 0) ? 3 : tetraD.x;
	}
	float2 uvRotA = RotateUV(uv, thetaBins[tetraA.x] * DEG2RAD, 0.0.rr);
	float2 uvRotB = RotateUV(uv, thetaBins[tetraB.x] * DEG2RAD, 0.0.rr);
	float2 uvRotC = RotateUV(uv, thetaBins[tetraC.x] * DEG2RAD, 0.0.rr);
	float2 uvRotD = RotateUV(uv, thetaBins[tetraD.x] * DEG2RAD, 0.0.rr);

	// SAMPLE GLINT GRIDS
	uint gridSeedA = asuint(HashWithoutSine13(float3(log2(divLods[tetraA.z]), thetaBins[tetraA.x] % 360, ratios[tetraA.y])));
	uint gridSeedB = asuint(HashWithoutSine13(float3(log2(divLods[tetraB.z]), thetaBins[tetraB.x] % 360, ratios[tetraB.y])));
	uint gridSeedC = asuint(HashWithoutSine13(float3(log2(divLods[tetraC.z]), thetaBins[tetraC.x] % 360, ratios[tetraC.y])));
	uint gridSeedD = asuint(HashWithoutSine13(float3(log2(divLods[tetraD.z]), thetaBins[tetraD.x] % 360, ratios[tetraD.y])));
	float4 sampleA = SampleGlintGridSimplexMultinomial(uvRotA / divLods[tetraA.z] / float2(1.0, ratios[tetraA.y]), gridSeedA, ratios[tetraA.y] * footprintAreas[tetraA.z], successProbs, tetraBarycentricWeights.x);
	float4 sampleB = SampleGlintGridSimplexMultinomial(uvRotB / divLods[tetraB.z] / float2(1.0, ratios[tetraB.y]), gridSeedB, ratios[tetraB.y] * footprintAreas[tetraB.z], successProbs, tetraBarycentricWeights.y);
	float4 sampleC = SampleGlintGridSimplexMultinomial(uvRotC / divLods[tetraC.z] / float2(1.0, ratios[tetraC.y]), gridSeedC, ratios[tetraC.y] * footprintAreas[tetraC.z], successProbs, tetraBarycentricWeights.z);
	float4 sampleD = SampleGlintGridSimplexMultinomial(uvRotD / divLods[tetraD.z] / float2(1.0, ratios[tetraD.y]), gridSeedD, ratios[tetraD.y] * footprintAreas[tetraD.z], successProbs, tetraBarycentricWeights.w);
	float4 result = sampleA + sampleB + sampleC + sampleD;
	result = SanitizePositiveFinite(result);
	return result;
}
