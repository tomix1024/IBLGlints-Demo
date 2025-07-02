
// Include Deliot & Belcour implementation
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/GlintsMultinomial.hlsl"

float2x4 SampleGlintGridSimplexVertexMultinomial8(float2x4 randSlopesB, float2x4 randSlopesG, float2x4 successProbs, float microfacetCountBlended)
{
	// NOTE: Require successProbs to be normalized!
	float successProbA = successProbs[0].x + successProbs[0].y + successProbs[0].z + successProbs[0].w;
	float successProbB = max(0, 1 - successProbA);
	float4 successProbXinA = successProbs[0] / max(1e-8, successProbA);
	float4 successProbXinB = successProbs[1] / max(1e-8, successProbB);

	float result_A, result_B;
	SampleBinomialValueForSurfaceCellMultinomial(randSlopesB[0].w, randSlopesG[0].w, successProbA, microfacetCountBlended, result_A, result_B);
	//result_A = successProbA * microfacetCountBlended;
	//result_B = successProbB * microfacetCountBlended;

	float4 result_0 = SampleGlintGridSimplexVertexMultinomial(randSlopesB[0], randSlopesG[0], successProbXinA, result_A);
	float4 result_1 = SampleGlintGridSimplexVertexMultinomial(randSlopesB[1], randSlopesG[1], successProbXinB, result_B);

	float2x4 result;
	result[0] = result_0;
	result[1] = result_1;
	return result;
}

float2x4 SampleGlintGridSimplexMultinomial8(float2 uv, uint gridSeed, float footprintArea, float2x4 successProbs, float gridWeight)
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
	float3 rand0x = UintToUniformFloat(seed0); // pcg3dFloat(uint3(glint0 + 2147483648, gridSeed)); // TODO : optimize away manual seeds
	float3 rand1x = UintToUniformFloat(seed1); // pcg3dFloat(uint3(glint1 + 2147483648, gridSeed));
	float3 rand2x = UintToUniformFloat(seed2); // pcg3dFloat(uint3(glint2 + 2147483648, gridSeed));
	seed0 = pcg3d(seed0);
	seed1 = pcg3d(seed1);
	seed2 = pcg3d(seed2);
	float3 rand0y = UintToUniformFloat(seed0); // pcg3dFloat(uint3(glint0 + 2147483648, gridSeed)); // TODO : optimize away manual seeds
	float3 rand1y = UintToUniformFloat(seed1); // pcg3dFloat(uint3(glint1 + 2147483648, gridSeed));
	float3 rand2y = UintToUniformFloat(seed2); // pcg3dFloat(uint3(glint2 + 2147483648, gridSeed));

	// Random numbers for footprintMicrofacetCount
	float3 rand0z = pcg3dFloat(seed0);
	float3 rand1z = pcg3dFloat(seed1);
	float3 rand2z = pcg3dFloat(seed2);

	// Compute microfacet count with randomization
	float3 logDensityRand = clamp(SampleNormalDistribution(float3(rand0x.x, rand1x.x, rand2x.x), _LogMicrofacetDensity, _DensityRandomization), 0.0, 50.0); // TODO : optimize sampleNormalDist
	float3 microfacetCountTotal = exp(logDensityRand);
	float3 microfacetCount = microfacetCountTotal * gridWeight;
	float3 microfacetCountBlended = ComputeMicrofacetCountInFootprint(rand0y, rand1y, rand2y, footprintArea, microfacetCount);
	float3 microfacetCountAverage = microfacetCount*footprintArea;

	// Generate per surface cell random numbers
	float2x4 rand0SlopesB, rand1SlopesB, rand2SlopesB;
	float2x4 rand0SlopesG, rand1SlopesG, rand2SlopesG;
	float2 slopeLerpDummy;
	CustomRand4Texture(float2(0,0), rand0x.yz, rand0SlopesB[0], rand0SlopesG[0], slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand1x.yz, rand1SlopesB[0], rand1SlopesG[0], slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand2x.yz, rand2SlopesB[0], rand2SlopesG[0], slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand0z.yz, rand0SlopesB[1], rand0SlopesG[1], slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand1z.yz, rand1SlopesB[1], rand1SlopesG[1], slopeLerpDummy);
	CustomRand4Texture(float2(0,0), rand2z.yz, rand2SlopesB[1], rand2SlopesG[1], slopeLerpDummy);

	// Generate numbers of reflecting microfacets
	float2x4 result0 = SampleGlintGridSimplexVertexMultinomial8(rand0SlopesB, rand0SlopesG, successProbs, microfacetCountBlended.x);
	float2x4 result1 = SampleGlintGridSimplexVertexMultinomial8(rand1SlopesB, rand1SlopesG, successProbs, microfacetCountBlended.y);
	float2x4 result2 = SampleGlintGridSimplexVertexMultinomial8(rand2SlopesB, rand2SlopesG, successProbs, microfacetCountBlended.z);

	// Divide by each microfacet count by expected microfacet count inside of grid cell
	// NOTE: NO GRID WEIGHT IN THIS DIVISION!
	result0 = result0 / (microfacetCountTotal.x*footprintArea);
	result1 = result1 / (microfacetCountTotal.y*footprintArea);
	result2 = result2 / (microfacetCountTotal.z*footprintArea);

	// Interpolate across surface cells
	float2x4 result = result0 * barycentrics.x + result1 * barycentrics.y + result2 * barycentrics.z;
	return result;
}

float2x4 SampleGlintsMultinomial8NDF_Internal(float2x4 successProbs, float2 uv, float2 duvdx, float2 duvdy)
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
	float2x4 sampleA = SampleGlintGridSimplexMultinomial8(uvRotA / divLods[tetraA.z] / float2(1.0, ratios[tetraA.y]), gridSeedA, ratios[tetraA.y] * footprintAreas[tetraA.z], successProbs, tetraBarycentricWeights.x);
	float2x4 sampleB = SampleGlintGridSimplexMultinomial8(uvRotB / divLods[tetraB.z] / float2(1.0, ratios[tetraB.y]), gridSeedB, ratios[tetraB.y] * footprintAreas[tetraB.z], successProbs, tetraBarycentricWeights.y);
	float2x4 sampleC = SampleGlintGridSimplexMultinomial8(uvRotC / divLods[tetraC.z] / float2(1.0, ratios[tetraC.y]), gridSeedC, ratios[tetraC.y] * footprintAreas[tetraC.z], successProbs, tetraBarycentricWeights.z);
	float2x4 sampleD = SampleGlintGridSimplexMultinomial8(uvRotD / divLods[tetraD.z] / float2(1.0, ratios[tetraD.y]), gridSeedD, ratios[tetraD.y] * footprintAreas[tetraD.z], successProbs, tetraBarycentricWeights.w);
	float2x4 result = sampleA + sampleB + sampleC + sampleD;
	result = SanitizePositiveFinite(result);
	return result;
}
