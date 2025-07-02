#ifndef GLINT_RANDOM_HLSL
#define GLINT_RANDOM_HLSL

#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/PCG3D.hlsl"
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/HashWithoutSine.hlsl"

float UintToUniformFloat(uint v)
{
    //float res = float(v) * (1.0 / 4294967296.0);
    // First generate a random number in [1, 2) and subtract 1.
    uint mantissa_bits = v >> 9;
    // For the paranoid user that wants to ensure random numbers in the open interval (1, 2):
    //uint mantissa_bits = (rngState & (~0 << 9)) != 0 ? rngState >> 9 : ((rngState << 14) | (1 << 13));
    uint float_bits = mantissa_bits | 0x3f800000;
    float res = asfloat(float_bits) - 1.0;
    return res;
}

float2 UintToUniformFloat(uint2 v)
{
    //float res = float(v) * (1.0 / 4294967296.0);
    // First generate a random number in [1, 2) and subtract 1.
    uint2 mantissa_bits = v >> 9;
    // For the paranoid user that wants to ensure random numbers in the open interval (1, 2):
    //uint mantissa_bits = (rngState & (~0 << 9)) != 0 ? rngState >> 9 : ((rngState << 14) | (1 << 13));
    uint2 float_bits = mantissa_bits | 0x3f800000;
    float2 res = float2(asfloat(float_bits.x), asfloat(float_bits.y)) - 1.0;
    return res;
}

float3 UintToUniformFloat(uint3 v)
{
    //float res = float(v) * (1.0 / 4294967296.0);
    // First generate a random number in [1, 2) and subtract 1.
    uint3 mantissa_bits = v >> 9;
    // For the paranoid user that wants to ensure random numbers in the open interval (1, 2):
    //uint mantissa_bits = (rngState & (~0 << 9)) != 0 ? rngState >> 9 : ((rngState << 14) | (1 << 13));
    uint3 float_bits = mantissa_bits | 0x3f800000;
    float3 res = asfloat(float_bits) - 1.0;
    return res;
}

float4 UintToUniformFloat(uint4 v)
{
    //float res = float(v) * (1.0 / 4294967296.0);
    // First generate a random number in [1, 2) and subtract 1.
    uint4 mantissa_bits = v >> 9;
    // For the paranoid user that wants to ensure random numbers in the open interval (1, 2):
    //uint mantissa_bits = (rngState & (~0 << 9)) != 0 ? rngState >> 9 : ((rngState << 14) | (1 << 13));
    uint4 float_bits = mantissa_bits | 0x3f800000;
    float4 res = float4(asfloat(float_bits.x), asfloat(float_bits.y), asfloat(float_bits.z), asfloat(float_bits.w)) - 1.0;
    return res;
}

float2 pcg2dFloat(uint2 v)
{
    v = pcg2d(v);
    return UintToUniformFloat(v);
}

float3 pcg3dFloat(uint3 v)
{
    v = pcg3d(v);
    return UintToUniformFloat(v);
}

float4 pcg4dFloat(uint4 v)
{
    v = pcg4d(v);
    return UintToUniformFloat(v);
}


float ErfInv(float x)
{
	float w, p;
	w = -log((1.0 - x) * (1.0 + x));
	if (w < 5.000000)
	{
		w = w - 2.500000;
		p = 2.81022636e-08;
		p = 3.43273939e-07 + p * w;
		p = -3.5233877e-06 + p * w;
		p = -4.39150654e-06 + p * w;
		p = 0.00021858087 + p * w;
		p = -0.00125372503 + p * w;
		p = -0.00417768164 + p * w;
		p = 0.246640727 + p * w;
		p = 1.50140941 + p * w;
	}
	else
	{
		w = sqrt(w) - 3.000000;
		p = -0.000200214257;
		p = 0.000100950558 + p * w;
		p = 0.00134934322 + p * w;
		p = -0.00367342844 + p * w;
		p = 0.00573950773 + p * w;
		p = -0.0076224613 + p * w;
		p = 0.00943887047 + p * w;
		p = 1.00167406 + p * w;
		p = 2.83297682 + p * w;
	}
	return p * x;
}

float SampleNormalDistribution(float u, float mu, float sigma)
{
	//return mu + sigma * (sqrt(-2.0 * log(u.x))* cos(2.0 * pi * u.y));
	float x0 = sigma * 1.414213f * ErfInv(2.0 * u - 1.0) + mu;
	return x0;
}

float3 SampleNormalDistribution(float3 u, float mu, float sigma)
{
	//return mu + sigma * (sqrt(-2.0 * log(u.x))* cos(2.0 * pi * u.y));
	float x0 = sigma * 1.414213f * ErfInv(2.0 * u.x - 1.0) + mu;
	float x1 = sigma * 1.414213f * ErfInv(2.0 * u.y - 1.0) + mu;
	float x2 = sigma * 1.414213f * ErfInv(2.0 * u.z - 1.0) + mu;
	return float3(x0, x1, x2);
}

float4 SampleNormalDistribution(float4 u, float mu, float sigma)
{
	//return mu + sigma * (sqrt(-2.0 * log(u.x))* cos(2.0 * pi * u.y));
	float x0 = sigma * 1.414213f * ErfInv(2.0 * u.x - 1.0) + mu;
	float x1 = sigma * 1.414213f * ErfInv(2.0 * u.y - 1.0) + mu;
	float x2 = sigma * 1.414213f * ErfInv(2.0 * u.z - 1.0) + mu;
	float x3 = sigma * 1.414213f * ErfInv(2.0 * u.w - 1.0) + mu;
	return float4(x0, x1, x2, x3);
}


#endif // GLINT_RANDOM_HLSL
