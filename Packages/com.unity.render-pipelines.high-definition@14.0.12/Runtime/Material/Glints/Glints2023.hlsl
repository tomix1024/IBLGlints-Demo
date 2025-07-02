#ifndef GLINTS2023_HLSL
#define GLINTS2023_HLSL

#define m_pi 3.141592
#define m_i_pi 0.318309
#define m_i_sqrt_2 0.707106
static const float DEG2RAD = 0.01745329251;
static const float RAD2DEG = 57.2957795131;

Texture2D<float4> _Glint2023NoiseMap;
float4 _Glint2023NoiseMap_TexelSize;
float _ScreenSpaceScale;
float _LogMicrofacetDensity;
float _HalfwaySlopeScale; // formerly == 1 / _MicrofacetRoughness
float _DensityRandomization;
// Use uniform variabels for debugging...
float _FixGatingNlt1;


float Remap(float s, float a1, float a2, float b1, float b2)
{
	return b1 + (s - a1) * (b2 - b1) / (a2 - a1);
}

float Remap01To(float s, float b1, float b2)
{
	return b1 + s * (b2 - b1);
}

float RemapTo01(float s, float a1, float a2)
{
	return (s - a1) / (a2 - a1);
}

float2 RemapTo01(float2 s, float2 a1, float2 a2)
{
	return (s - a1) / (a2 - a1);
}

float4 RemapTo01(float4 s, float4 a1, float4 a2)
{
	return (s - a1) / (a2 - a1);
}

#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Common.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Random.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Binomial.hlsl"


//=======================================================================================
// TOOLS
//=======================================================================================

float2x2 Inverse(float2x2 A)
{
	return float2x2(A[1][1], -A[0][1], -A[1][0], A[0][0]) / max(1e-6, determinant(A));
}

void GetGradientEllipse(float2 duvdx, float2 duvdy, out float2 ellipseMajorAxis, out float2 ellipseMinorAxis, out float ellipseMajorLength, out float ellipseMinorLength)
{
	float2x2 J = float2x2(duvdx, duvdy);
	J = Inverse(J);
	J = mul(J, transpose(J));

	float a = J[0][0];
	float b = J[0][1];
	float c = J[1][0];
	float d = J[1][1];

	float T = a + d;
	float D = a * d - b * c;
	float L1 = T / 2.0 - sqrt(max(0, T * T / 3.99999 - D));
	float L2 = T / 2.0 + sqrt(max(0, T * T / 3.99999 - D));

	float2 A0 = float2(L1 - d, c);
	float2 A1 = float2(L2 - d, c);
	float r0 = rsqrt(L1);
	float r1 = rsqrt(L2);
	//ellipseMajor = normalize(A0) * r0;
	//ellipseMinor = normalize(A1) * r1;
	ellipseMajorAxis = normalize(A0);
	ellipseMinorAxis = normalize(A1);
	ellipseMajorLength = r0;
	ellipseMinorLength = r1;
}

void GetGradientEllipseNew(float2 duvdx, float2 duvdy, out float2 ellipseMajorAxis, out float2 ellipseMinorAxis, out float ellipseMajorLength, out float ellipseMinorLength)
{
	float2x2 J = float2x2(duvdx, duvdy);
	// Do computation without inversion! (=> eigen values are inverted)
	//J = Inverse(J);
	// Also communte J.t @ J ...
	J = mul(transpose(J), J);

	// NOTE: J is symmetric here, therefore b=c
	float a = J[0][0];
	float b = J[0][1];
	float c = b;//J[1][0];
	float d = J[1][1];

	float T = a + d;
	float D = a * d - b * c;
	float L1 = T / 2.0 - sqrt(max(0, T * T / 4 - D));
	float L2 = T / 2.0 + sqrt(max(0, T * T / 4 - D));

	float2 A0 = float2(L1 - d, c);//lerp(float2(L1 - d, c), float2(b, L1 - a), abs(L1 - a) > abs(L1 - d));
	float2 A1 = float2(L2 - d, c);//lerp(float2(L2 - d, c), float2(b, L2 - a), abs(L2 - a) > abs(L2 - d));
	float r0 = sqrt(L1);
	float r1 = sqrt(L2);
	//ellipseMinor = normalize(A0) * r0;
	//ellipseMajor = normalize(A1) * r1;
	ellipseMajorAxis = -normalize(A1);
	ellipseMinorAxis = -normalize(A0);
	ellipseMajorLength = r1;
	ellipseMinorLength = r0;
}

float2 RotateUV(float2 uv, float rotation, float2 mid)
{
	return float2(
		cos(rotation) * (uv.x - mid.x) + sin(rotation) * (uv.y - mid.y) + mid.x,
		cos(rotation) * (uv.y - mid.y) - sin(rotation) * (uv.x - mid.x) + mid.y
		);
}

float BilinearLerp(float4 values, float2 valuesLerp)
{
	// Values XY = float4(00, 01, 10, 11)
	float resultX = lerp(values.x, values.z, valuesLerp.x);
	float resultY = lerp(values.y, values.w, valuesLerp.x);
	float result = lerp(resultX, resultY, valuesLerp.y);
	return result;
}

float4 BilinearLerpParallel4(float4 values00, float4 values01, float4 values10, float4 values11, float4 valuesLerpX, float4 valuesLerpY)
{
	// Values XY = float4(00, 01, 10, 11)
	float4 resultX = lerp(values00, values10, valuesLerpX);
	float4 resultY = lerp(values01, values11, valuesLerpX);
	float4 result = lerp(resultX, resultY, valuesLerpY);
	return result;
}

float3 GetBarycentricWeights(float2 p, float2 a, float2 b, float2 c)
{
	/*float2 v0 = b - a;
	float2 v1 = c - a;
	float2 v2 = p - a;
	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d11 = dot(v1, v1);
	float d20 = dot(v2, v0);
	float d21 = dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0 - v - w;
	return float3(u, v, w);*/

	float2 v0 = b - a;
	float2 v1 = c - a;
	float2 v2 = p - a;
	float den = v0.x * v1.y - v1.x * v0.y;
	float v = (v2.x * v1.y - v1.x * v2.y) / den;
	float w = (v0.x * v2.y - v2.x * v0.y) / den;
	float u = 1.0f - v - w;
	return float3(u, v, w);
}

float4 GetBarycentricWeightsTetrahedron(float3 p, float3 v1, float3 v2, float3 v3, float3 v4)
{
	float3 c11 = v1 - v4, c21 = v2 - v4, c31 = v3 - v4, c41 = v4 - p;

	float2 m1 = c31.yz / c31.x;
	float2 c12 = c11.yz - c11.x * m1, c22 = c21.yz - c21.x * m1, c32 = c41.yz - c41.x * m1;

	float4 uvwk = 0.0.rrrr;
	float m2 = c22.y / c22.x;
	uvwk.x = (c32.x * m2 - c32.y) / (c12.y - c12.x * m2);
	uvwk.y = -(c32.x + c12.x * uvwk.x) / c22.x;
	uvwk.z = -(c41.x + c21.x * uvwk.y + c11.x * uvwk.x) / c31.x;
	uvwk.w = 1.0 - uvwk.z - uvwk.y - uvwk.x;

	return uvwk;
}

void UnpackFloat(float input, out float a, out float b)
{
	// a = uniform uint16
	// b = gauss f16
	uint uintInput = asuint(input);
	a = (float(uintInput >> 16) + 0.5) / float(0xffff);
	b = f16tof32(uintInput);
}

void UnpackFloatParallel4(float4 input, out float4 a, out float4 b)
{
	// a = uniform uint16
	// b = gauss f16
	uint4 uintInput = asuint(input);
	a = (float4(uintInput >> 16) + 0.5) / float(0x10000);
	b = f16tof32(uintInput);
}


//=======================================================================================
// GLINTS TEST NOVEMBER 2022
//=======================================================================================
void CustomRand4Texture(float2 slope, float2 slopeRandOffset, out float4 outUniform, out float4 outGaussian, out float2 slopeLerp, int2 ioffset=int2(0,0))
{
	int2 size = _Glint2023NoiseMap_TexelSize.zw;
	float2 slope2 = abs(slope) * _HalfwaySlopeScale;
	slope2 = slope2 + (slopeRandOffset * size);
	slopeLerp = frac(slope2);
	int2 slopeCoord = (int2(floor(slope2)) + ioffset) % size;

	float4 packedRead = _Glint2023NoiseMap[slopeCoord];
	UnpackFloatParallel4(packedRead, outUniform, outGaussian);
}


float GenerateAngularBinomialValueForSurfaceCell(float4 randB, float4 randG, float2 slopeLerp, float footprintOneHitProba, float binomialSmoothWidth, float footprintMean, float footprintSTD, float microfacetCount)
{
	float4 gating;
	if (binomialSmoothWidth > 0.0000001)
		gating = saturate(RemapTo01(randB, footprintOneHitProba + binomialSmoothWidth, footprintOneHitProba - binomialSmoothWidth));
	else
		gating = randB < footprintOneHitProba;

	float4 gauss = randG * footprintSTD + footprintMean;
	gauss = clamp(gauss, 1, max(1, microfacetCount));
	float4 results = gating * gauss;
    float result = BilinearLerp(results, slopeLerp);
	return result;
}

void EvalSkewedGrid(float2 uv, out int2 vertex0, out int2 vertex1, out int2 vertex2, out float3 barycentrics)
{
	// Get surface space glint simplex grid cell
	const float2x2 gridToSkewedGrid = float2x2(1.0, -0.57735027, 0.0, 1.15470054);
	float2 skewedCoord = mul(gridToSkewedGrid, uv);
	int2 baseId = int2(floor(skewedCoord));
	float3 temp = float3(frac(skewedCoord), 0.0);
	temp.z = 1.0 - temp.x - temp.y;
	float s = step(0.0, -temp.z);
	float s2 = 2.0 * s - 1.0;
	vertex0 = baseId + int2(s, s);
	vertex1 = baseId + int2(s, 1.0 - s);
	vertex2 = baseId + int2(1.0 - s, s);
	barycentrics = float3(-temp.z * s2, s - temp.y * s2, s - temp.x * s2);
}


float3 ComputeMicrofacetCountInFootprint(float3 rand0, float3 rand1, float3 rand2, float footprintArea, float3 microfacetCount)
{
	float3 randU0 = float3(rand0.x, rand1.x, rand2.x);
	float3 randU1 = float3(rand0.y, rand1.y, rand2.y);
	float3 randU2 = float3(rand0.z, rand1.z, rand2.z);
	// TODO eval box-mueller vs erfinv approx speed!
	float3 randG1 = SampleNormalDistribution(randU1, 0, 1);

	// NOTE: Assume footprintArea < 1

	float3 footprintMicrofacetCount;
#if defined(_GLINTS_SURFACE_DISTRIBUTION_DUALGATED)
	{
		// Dual-gated Gaussian
		float4 NfPos;
		float4 NfNeg;
		// TODO implement for float3...
		// In the compiler we trust!
		EvaluateBinomialValueDualGated(randU0.xyzz, randG1.xyzz, footprintArea, microfacetCount.xyzz, NfPos, NfNeg);
		footprintMicrofacetCount = NfPos.xyz;
	}
#elif defined(_GLINTS_SURFACE_DISTRIBUTION_GATED)
	{
		// Single-gated Gaussian
		footprintMicrofacetCount = EvaluateBinomialValueGated(randU0.xyzz, randG1.xyzz, footprintArea, microfacetCount.xyzz).xyz;
	}
#elif defined(_GLINTS_SURFACE_DISTRIBUTION_UNGATED)
	{
		// Ungated Gaussian
		footprintMicrofacetCount = EvaluateBinomialValueUngated(randG1.xyzz, footprintArea, microfacetCount.xyzz).xyz;
	}
#else
	{
		// Uniform
		footprintMicrofacetCount = microfacetCount * footprintArea;
	}
#endif
	// Rounding happens implicitly during evaluation!
	return footprintMicrofacetCount;
}



float SampleGlintGridSimplex(float2 uv, uint gridSeed, float2 slope, float footprintArea, float targetNDF, float gridWeight)
{
	int2 glint0, glint1, glint2;
	float3 barycentrics;
	EvalSkewedGrid(uv, glint0, glint1, glint2, barycentrics);

	// Generate per surface cell random numbers
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

	// Get per surface cell per slope cell random numbers
	float4 rand0SlopesB, rand1SlopesB, rand2SlopesB;
	float4 rand0SlopesG, rand1SlopesG, rand2SlopesG;
	float2 slopeLerp0, slopeLerp1, slopeLerp2;
	CustomRand4Texture(slope, rand0.yz, rand0SlopesB, rand0SlopesG, slopeLerp0);
	CustomRand4Texture(slope, rand1.yz, rand1SlopesB, rand1SlopesG, slopeLerp1);
	CustomRand4Texture(slope, rand2.yz, rand2SlopesB, rand2SlopesG, slopeLerp2);

	// Compute binomial properties
	float hitProba = targetNDF; // probability of hitting desired half vector in NDF distribution
	float3 footprintOneHitProba = (1.0 - OneMinusP_Pow_N(hitProba.rrr, microfacetCountBlended)); // probability of hitting at least one microfacet in footprint
	float3 footprintMean = 1 + max(0, microfacetCountBlended - 1.0) * hitProba.rrr; // Expected value of number of hits in the footprint given already one hit
	float3 footprintSTD = sqrt(max(0, microfacetCountBlended - 1.0) * hitProba.rrr * (1.0 - hitProba.rrr)); // Standard deviation of number of hits in the footprint given already one hit

	// This does the trick! It's even continuous around N=1
	footprintOneHitProba = lerp(footprintOneHitProba, hitProba * microfacetCountBlended, (microfacetCountBlended < 1) * _FixGatingNlt1);

	float3 binomialSmoothWidth = 0.1 * clamp(footprintOneHitProba * 10, 0.0, 1.0) * clamp((1.0 - footprintOneHitProba) * 10, 0.0, 1.0);

	// Generate numbers of reflecting microfacets
	float result0, result1, result2;
	result0 = GenerateAngularBinomialValueForSurfaceCell(rand0SlopesB, rand0SlopesG, slopeLerp0, footprintOneHitProba.x, binomialSmoothWidth.x, footprintMean.x, footprintSTD.x, microfacetCountBlended.x);
	result1 = GenerateAngularBinomialValueForSurfaceCell(rand1SlopesB, rand1SlopesG, slopeLerp1, footprintOneHitProba.y, binomialSmoothWidth.y, footprintMean.y, footprintSTD.y, microfacetCountBlended.y);
	result2 = GenerateAngularBinomialValueForSurfaceCell(rand2SlopesB, rand2SlopesG, slopeLerp2, footprintOneHitProba.z, binomialSmoothWidth.z, footprintMean.z, footprintSTD.z, microfacetCountBlended.z);

	// Interpolate result for glint grid cell
	// TODO make division fix optional via GUI! (same for Multinomial)
	float3 results = float3(result0, result1, result2) / (microfacetCountTotal*footprintArea);
	float result = dot(results, barycentrics);
	return result;
}

void GetAnisoCorrectingGridTetrahedron(bool centerSpecialCase, inout float thetaBinLerp, float ratioLerp, float lodLerp, out float3 p0, out float3 p1, out float3 p2, out float3 p3)
{
	[branch] if (centerSpecialCase == true) // SPECIAL CASE (no anisotropy, center of blending pattern, different triangulation)
	{
		float3 a = float3(0, 1, 0);
		float3 b = float3(0, 0, 0);
		float3 c = float3(1, 1, 0);
		float3 d = float3(0, 1, 1);
		float3 e = float3(0, 0, 1);
		float3 f = float3(1, 1, 1);
		[branch] if (lodLerp > 1.0 - ratioLerp) // Upper pyramid
		{
			[branch] if (RemapTo01(lodLerp, 1.0 - ratioLerp, 1.0) > thetaBinLerp) // Left-up tetrahedron (a, e, d, f)
			{
				p0 = a; p1 = e; p2 = d; p3 = f;
			}
			else // Right-down tetrahedron (f, e, c, a)
			{
				p0 = f; p1 = e; p2 = c; p3 = a;
			}
		}
		else // Lower tetrahedron (b, a, c, e)
		{
			p0 = b; p1 = a; p2 = c; p3 = e;
		}
	}
	else // NORMAL CASE
	{
		float3 a = float3(0, 1, 0);
		float3 b = float3(0, 0, 0);
		float3 c = float3(0.5, 1, 0);
		float3 d = float3(1, 0, 0);
		float3 e = float3(1, 1, 0);
		float3 f = float3(0, 1, 1);
		float3 g = float3(0, 0, 1);
		float3 h = float3(0.5, 1, 1);
		float3 i = float3(1, 0, 1);
		float3 j = float3(1, 1, 1);
		[branch] if (thetaBinLerp < 0.5 && thetaBinLerp * 2.0 < ratioLerp) // Prism A
		{
			[branch] if (lodLerp > 1.0 - ratioLerp) // Upper pyramid
			{
				[branch] if (RemapTo01(lodLerp, 1.0 - ratioLerp, 1.0) > RemapTo01(thetaBinLerp * 2.0, 0.0, ratioLerp)) // Left-up tetrahedron (a, f, h, g)
				{
					p0 = a; p1 = f; p2 = h; p3 = g;
				}
				else // Right-down tetrahedron (c, a, h, g)
				{
					p0 = c; p1 = a; p2 = h; p3 = g;
				}
			}
			else // Lower tetrahedron (b, a, c, g)
			{
				p0 = b; p1 = a; p2 = c; p3 = g;
			}
		}
		else if (1.0 - ((thetaBinLerp - 0.5) * 2.0) > ratioLerp) // Prism B
		{
			[branch] if (lodLerp < 1.0 - ratioLerp) // Lower pyramid
			{
				[branch] if (RemapTo01(lodLerp, 0.0, 1.0 - ratioLerp) > RemapTo01(thetaBinLerp, 0.5 - (1.0 - ratioLerp) * 0.5, 0.5 + (1.0 - ratioLerp) * 0.5)) // Left-up tetrahedron (b, g, i, c)
				{
					p0 = b; p1 = g; p2 = i; p3 = c;
				}
				else // Right-down tetrahedron (d, b, c, i)
				{
					p0 = d; p1 = b; p2 = c; p3 = i;
				}
			}
			else // Upper tetrahedron (c, g, h, i)
			{
				p0 = c; p1 = g; p2 = h; p3 = i;
			}
		}
		else // Prism C
		{
			[branch] if (lodLerp > 1.0 - ratioLerp) // Upper pyramid
			{
				[branch] if (RemapTo01(lodLerp, 1.0 - ratioLerp, 1.0) > RemapTo01((thetaBinLerp - 0.5) * 2.0, 1.0 - ratioLerp, 1.0)) // Left-up tetrahedron (c, j, h, i)
				{
					p0 = c; p1 = j; p2 = h; p3 = i;
				}
				else // Right-down tetrahedron (e, i, c, j)
				{
					p0 = e; p1 = i; p2 = c; p3 = j;
				}
			}
			else // Lower tetrahedron (d, e, c, i)
			{
				p0 = d; p1 = e; p2 = c; p3 = i;
			}
		}
	}

	return;
}

float SampleGlints2023NDF_Internal(float3 localHalfVector, float successProb, float2 uv, float2 duvdx, float2 duvdy)
{
	// ACCURATE PIXEL FOOTPRINT ELLIPSE
	float2 ellipseMajorAxis, ellipseMinorAxis;
	float ellipseMajorLength, ellipseMinorLength;
	GetGradientEllipseNew(duvdx, duvdy, ellipseMajorAxis, ellipseMinorAxis, ellipseMajorLength, ellipseMinorLength);
	float ellipseRatio = ellipseMajorLength / ellipseMinorLength;

	// SHARED GLINT NDF VALUES
	float2 slope = localHalfVector.xy; // Orthogrtaphic slope projected grid

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
	float sampleA = SampleGlintGridSimplex(uvRotA / divLods[tetraA.z] / float2(1.0, ratios[tetraA.y]), gridSeedA, slope, ratios[tetraA.y] * footprintAreas[tetraA.z], successProb, tetraBarycentricWeights.x);
	float sampleB = SampleGlintGridSimplex(uvRotB / divLods[tetraB.z] / float2(1.0, ratios[tetraB.y]), gridSeedB, slope, ratios[tetraB.y] * footprintAreas[tetraB.z], successProb, tetraBarycentricWeights.y);
	float sampleC = SampleGlintGridSimplex(uvRotC / divLods[tetraC.z] / float2(1.0, ratios[tetraC.y]), gridSeedC, slope, ratios[tetraC.y] * footprintAreas[tetraC.z], successProb, tetraBarycentricWeights.z);
	float sampleD = SampleGlintGridSimplex(uvRotD / divLods[tetraD.z] / float2(1.0, ratios[tetraD.y]), gridSeedD, slope, ratios[tetraD.y] * footprintAreas[tetraD.z], successProb, tetraBarycentricWeights.w);
	return SanitizePositiveFinite(sampleA + sampleB + sampleC + sampleD);
}

#endif // GLINTS2023_HLSL
