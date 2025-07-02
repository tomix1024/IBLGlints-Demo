#ifndef GLINT_BINOMIAL_HLSL
#define GLINT_BINOMIAL_HLSL

// Global uniform parameter to control "smoothness" of gating transition
float _GatingSmoothWidth;

float EvalGating(float rand, float prob)
{
	float smoothWidth = _GatingSmoothWidth * saturate(prob / _GatingSmoothWidth) * saturate((1.0 - prob) / _GatingSmoothWidth);
	float gating;
	if (smoothWidth > 0.0000001)
		gating = saturate(RemapTo01(rand, prob + smoothWidth, prob - smoothWidth));
	else
		gating = rand < prob;
	return gating;
}

float2 EvalGating(float2 rand, float prob)
{
	float smoothWidth = _GatingSmoothWidth * saturate(prob / _GatingSmoothWidth) * saturate((1.0 - prob) / _GatingSmoothWidth);
	float2 gating;
	if (smoothWidth > 0.0000001)
		gating = saturate(RemapTo01(rand, prob + smoothWidth, prob - smoothWidth));
	else
		gating = rand < prob;
	return gating;
}

float2 EvalGating(float2 rand, float2 prob)
{
	float2 smoothWidth = _GatingSmoothWidth * saturate(prob / _GatingSmoothWidth) * saturate((1.0 - prob) / _GatingSmoothWidth);
	float2 gatingSmooth = saturate(RemapTo01(rand, prob + smoothWidth, prob - smoothWidth));
	float2 gatingHard = rand < prob;
	float2 gating = lerp(gatingHard, gatingSmooth, smoothWidth > 0.0000001);
	return gating;
}

float4 EvalGating(float4 rand, float prob)
{
	float smoothWidth = _GatingSmoothWidth * saturate(prob / _GatingSmoothWidth) * saturate((1.0 - prob) / _GatingSmoothWidth);
	float4 gating;
	if (smoothWidth > 0.0000001)
		gating = saturate(RemapTo01(rand, prob + smoothWidth, prob - smoothWidth));
	else
		gating = rand < prob;
	return gating;
}

float4 EvalGating(float4 rand, float4 prob)
{
	float4 smoothWidth = _GatingSmoothWidth * saturate(prob / _GatingSmoothWidth) * saturate((1.0 - prob) / _GatingSmoothWidth);
	float4 gatingSmooth = saturate(RemapTo01(rand, prob + smoothWidth, prob - smoothWidth));
	float4 gatingHard = rand < prob;
	float4 gating = lerp(gatingHard, gatingSmooth, smoothWidth > 0.0000001);
	return gating;
}

//// OneMinusPPowN

float OneMinusP_Pow_N(float p, float N, float eps=0.00028840315031266055)
{
	// Problem: (1-p) rounds up to 1 for small p.
	// This becomes very noticeable for p < 1e-7, and measurable even for p = 5e-4...
	// For large N, (1-p)^N might still yield some value in (0, 1).
	// Workaround:
	// Consider [ (1-p)^2 ]^(N/2) = [1 + p^2 - 2p]^(N/2)
	// Notice that p^2 is even smaller than p.
	// We replace 2 with some constant c here and assume p^2 = 0.
	// For epsilon = 1e-3.54 we can accet (1-epsilon)^(N/c)
	// This requires epsilon = c*p, i.e. c = epsilon/p.
	// For p > epsilon we simply set c=1, such that the original formula is not modified!
	// Note that p might actually approach zero such that c is unbounded, which also leads to errors!
	float c = clamp(eps/p, 1, 1e6);
#ifdef _GLINTS_BROKEN_ONEMINUSPPOWN
	c = 1;
#endif // _GLINTS_BROKEN_ONEMINUSPPOWN
	float result = pow(1 - c*p, N/c);
	return result;
}

float2 OneMinusP_Pow_N(float2 p, float2 N, float eps=0.00028840315031266055)
{
	// Problem: (1-p) rounds up to 1 for small p.
	// This becomes very noticeable for p < 1e-7, and measurable even for p = 5e-4...
	// For large N, (1-p)^N might still yield some value in (0, 1).
	// Workaround:
	// Consider [ (1-p)^2 ]^(N/2) = [1 + p^2 - 2p]^(N/2)
	// Notice that p^2 is even smaller than p.
	// We replace 2 with some constant c here and assume p^2 = 0.
	// For epsilon = 1e-3.54 we can accet (1-epsilon)^(N/c)
	// This requires epsilon = c*p, i.e. c = epsilon/p.
	// For p > epsilon we simply set c=1, such that the original formula is not modified!
	// Note that p might actually approach zero such that c is unbounded, which also leads to errors!
	float2 c = clamp(eps/p, 1, 1e6);
#ifdef _GLINTS_BROKEN_ONEMINUSPPOWN
	c = 1.0.rr;
#endif // _GLINTS_BROKEN_ONEMINUSPPOWN
	float2 result = pow(1 - c*p, N/c);
	return result;
}

float3 OneMinusP_Pow_N(float3 p, float3 N, float eps=0.00028840315031266055)
{
	// Problem: (1-p) rounds up to 1 for small p.
	// This becomes very noticeable for p < 1e-7, and measurable even for p = 5e-4...
	// For large N, (1-p)^N might still yield some value in (0, 1).
	// Workaround:
	// Consider [ (1-p)^2 ]^(N/2) = [1 + p^2 - 2p]^(N/2)
	// Notice that p^2 is even smaller than p.
	// We replace 2 with some constant c here and assume p^2 = 0.
	// For epsilon = 1e-3.54 we can accet (1-epsilon)^(N/c)
	// This requires epsilon = c*p, i.e. c = epsilon/p.
	// For p > epsilon we simply set c=1, such that the original formula is not modified!
	// Note that p might actually approach zero such that c is unbounded, which also leads to errors!
	float3 c = clamp(eps/p, 1, 1e6);
#ifdef _GLINTS_BROKEN_ONEMINUSPPOWN
	c = 1.0.rrr;
#endif // _GLINTS_BROKEN_ONEMINUSPPOWN
	float3 result = pow(1 - c*p, N/c);
	return result;
}

float4 OneMinusP_Pow_N(float4 p, float4 N, float eps=0.00028840315031266055)
{
	// Problem: (1-p) rounds up to 1 for small p.
	// This becomes very noticeable for p < 1e-7, and measurable even for p = 5e-4...
	// For large N, (1-p)^N might still yield some value in (0, 1).
	// Workaround:
	// Consider [ (1-p)^2 ]^(N/2) = [1 + p^2 - 2p]^(N/2)
	// Notice that p^2 is even smaller than p.
	// We replace 2 with some constant c here and assume p^2 = 0.
	// For epsilon = 1e-3.54 we can accet (1-epsilon)^(N/c)
	// This requires epsilon = c*p, i.e. c = epsilon/p.
	// For p > epsilon we simply set c=1, such that the original formula is not modified!
	// Note that p might actually approach zero such that c is unbounded, which also leads to errors!
	float4 c = clamp(eps/p, 1, 1e6);
#ifdef _GLINTS_BROKEN_ONEMINUSPPOWN
	c = 1.0.rrrr;
#endif // _GLINTS_BROKEN_ONEMINUSPPOWN
	float4 result = pow(1 - c*p, N/c);
	return result;
}

//// EvaluateBinomialValueUngated

float EvaluateBinomialValueUngated(float randG, float successProb, float microfacetCount)
{
	float microfacetCountMean = microfacetCount * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = microfacetCount * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	float gauss = randG * microfacetCountStdDev + microfacetCountMean;
	float result = clamp(gauss, 0, microfacetCount);
	return result;
}

float2 EvaluateBinomialValueUngated(float2 randG, float successProb, float microfacetCount)
{
	float microfacetCountMean = microfacetCount * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = microfacetCount * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	float2 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	float2 result = clamp(gauss, 0, microfacetCount);
	return result;
}

float2 EvaluateBinomialValueUngated(float2 randG, float2 successProb, float2 microfacetCount)
{
	float2 microfacetCountMean = microfacetCount * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float2 microfacetCountVar = microfacetCount * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float2 microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	float2 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	float2 result = clamp(gauss, 0, microfacetCount);
	return result;
}

float4 EvaluateBinomialValueUngated(float4 randG, float successProb, float microfacetCount)
{
	float microfacetCountMean = microfacetCount * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = microfacetCount * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	float4 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	float4 result = clamp(gauss, 0, microfacetCount);
	return result;
}

float4 EvaluateBinomialValueUngated(float4 randG, float4 successProb, float4 microfacetCount)
{
	float4 microfacetCountMean = microfacetCount * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float4 microfacetCountVar = microfacetCount * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float4 microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	float4 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	float4 result = clamp(gauss, 0, microfacetCount);
	return result;
}

//// EvaluateBinomialValueGated

float EvaluateBinomialValueGated(float randB, float randG, float successProb, float microfacetCount)
{
	float microfacetCountOne = max(microfacetCount, 1);
	// Compute binomial properties
	float probOneHit = (1.0 - OneMinusP_Pow_N(successProb, microfacetCount)); // probability of hitting at least one microfacet in footprint
	float microfacetCountMean = 1 + (microfacetCountOne - 1.0) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = (microfacetCountOne - 1.0) * successProb * (1.0 - successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));

	// This does the trick! It's even continuous around N=1
	probOneHit = lerp(probOneHit, successProb * microfacetCount, (microfacetCount < 1) * _FixGatingNlt1);

	float gating = EvalGating(randB, probOneHit);

	float gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, microfacetCountOne);
	float result = gating * gauss;
	return result;
}

float2 EvaluateBinomialValueGated(float2 randB, float2 randG, float successProb, float microfacetCount)
{
	float microfacetCountOne = max(microfacetCount, 1);
	// Compute binomial properties
	float probOneHit = (1.0 - OneMinusP_Pow_N(successProb, microfacetCount)); // probability of hitting at least one microfacet in footprint
	float microfacetCountMean = 1 + (microfacetCountOne - 1.0) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = (microfacetCountOne - 1.0) * successProb * (1.0 - successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));

	// This does the trick! It's even continuous around N=1
	probOneHit = lerp(probOneHit, successProb * microfacetCount, (microfacetCount < 1) * _FixGatingNlt1);

	float2 gating = EvalGating(randB, probOneHit);

	float2 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, microfacetCountOne);
	float2 result = gating * gauss;
	return result;
}

float2 EvaluateBinomialValueGated(float2 randB, float2 randG, float2 successProb, float2 microfacetCount)
{
	float2 microfacetCountOne = max(microfacetCount, 1);
	// Compute binomial properties
	float2 probOneHit = (1.0 - OneMinusP_Pow_N(successProb, microfacetCount)); // probability of hitting at least one microfacet in footprint
	float2 microfacetCountMean = 1 + (microfacetCountOne - 1.0) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float2 microfacetCountVar = (microfacetCountOne - 1.0) * successProb * (1.0 - successProb); // Standard deviation of number of hits in the footprint given already one hit
	float2 microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));

	// This does the trick for microfacetCount < 1! It's even continuous around N=1
	probOneHit = lerp(probOneHit, successProb * microfacetCount, (microfacetCount < 1) * _FixGatingNlt1);

	float2 gating = EvalGating(randB, probOneHit);

	float2 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, microfacetCountOne);
	float2 result = gating * gauss;
	return result;
}

float4 EvaluateBinomialValueGated(float4 randB, float4 randG, float successProb, float microfacetCount)
{
	float microfacetCountOne = max(microfacetCount, 1);
	// Compute binomial properties
	float probOneHit = (1.0 - OneMinusP_Pow_N(successProb, microfacetCount)); // probability of hitting at least one microfacet in footprint
	float microfacetCountMean = 1 + (microfacetCountOne - 1.0) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = (microfacetCountOne - 1.0) * successProb * (1.0 - successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));

	// This does the trick! It's even continuous around N=1
	probOneHit = lerp(probOneHit, successProb * microfacetCount, (microfacetCount < 1) * _FixGatingNlt1);

	float4 gating = EvalGating(randB, probOneHit);

	float4 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, microfacetCountOne);
	float4 result = gating * gauss;
	return result;
}

float4 EvaluateBinomialValueGated(float4 randB, float4 randG, float4 successProb, float4 microfacetCount)
{
	float4 microfacetCountOne = max(microfacetCount, 1);
	// Compute binomial properties
	float4 probOneHit = (1.0 - OneMinusP_Pow_N(successProb, microfacetCount)); // probability of hitting at least one microfacet in footprint
	float4 microfacetCountMean = 1 + (microfacetCountOne - 1.0) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float4 microfacetCountVar = (microfacetCountOne - 1.0) * successProb * (1.0 - successProb); // Standard deviation of number of hits in the footprint given already one hit
	float4 microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));

	// This does the trick for microfacetCount < 1! It's even continuous around N=1
	probOneHit = lerp(probOneHit, successProb * microfacetCount, (microfacetCount < 1) * _FixGatingNlt1);

	float4 gating = EvalGating(randB, probOneHit);

	float4 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, microfacetCountOne);
	float4 result = gating * gauss;
	return result;
}

//// EvaluateBinomialValueDualGated

void EvaluateBinomialValueDualGated(float randB, float randG, float successProb, float microfacetCount, out float resultPos, out float resultNeg)
{
	// Unvectorized implementation

	successProb = saturate(successProb);
#ifdef ALLOW_FLIP
	bool flip = successProb > 0.5;
	if (flip) successProb = 1 - successProb;
#endif // ALLOW_FLIP

	// Evaluate probabilities for gating
	float probAllHitPos = saturate(microfacetCount - 1) * pow(successProb, max(2, microfacetCount));
	float probAllHitNeg = saturate(microfacetCount - 1) * OneMinusP_Pow_N(successProb, max(2, microfacetCount));
	float probExactOneHitPos = successProb * max(0, 1-abs(1-microfacetCount));
	float probExactOneHitNeg = (1-successProb) * max(0, 1-abs(1-microfacetCount));
	float probAllOrExactOneHitPos = probAllHitPos + probExactOneHitPos;
	float probAllOrExactOneHitNeg = probAllHitNeg + probExactOneHitNeg;

	// Evaluate gatings
	float gatingAllHitPos = EvalGating(randB, probAllHitPos);                   // randB < probAllHitPos;
	float gatingAllOrOneHitPos = EvalGating(randB, probAllOrExactOneHitPos);    // randB < probAllOrExactOneHitPos
	float gatingAllHitNeg = EvalGating(-randB, probAllHitNeg-1);                // randB > (1-probAllHitNeg);
	float gatingAllOrOneHitNeg = EvalGating(-randB, probAllOrExactOneHitNeg-1); // randB > (1-probAllOrExactOneHitNeg);
	//float gatingCentral = (microfacetCount > 1) && !(gatingAllOrOneHitPos || gatingAllOrOneHitNeg);
	float gatingCentral = (microfacetCount > 1) * (1-(gatingAllOrOneHitPos + gatingAllOrOneHitNeg));

	// Evaluate gaussian approx.
	float microfacetCountMean = 1 + (microfacetCount-2) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float microfacetCountVar = (microfacetCount-2) * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	// Transform to gaussian with mean + stddev.
	float gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, max(1, microfacetCount-1));

	//float gatingOneHitPos = gatingAllOrOneHitPos - gatingAllHitPos;
	//float gatingOneHitNeg = gatingAllOrOneHitNeg - gatingAllHitNeg;
	// resultPos = gatingOneHitPos * 1 + gatingAllHitPos * max(2, microfacetCount) + gatingCentral * gauss;
	// resultNeg = gatingOneHitNeg * 1 + gatingAllHitNeg * max(2, microfacetCount) + gatingCentral * (microfacetCount-gauss);
	resultPos = gatingAllOrOneHitPos + gatingAllHitPos * max(1, microfacetCount-1) + gatingCentral * gauss;
	resultNeg = gatingAllOrOneHitNeg + gatingAllHitNeg * max(1, microfacetCount-1) + gatingCentral * (max(2, microfacetCount)-gauss);
	//resultPos = lerp(gatingAllOrOneHitPos, max(2, microfacetCount), gatingAllHitPos) + gatingCentral * gauss;
	//resultNeg = lerp(gatingAllOrOneHitNeg, max(2, microfacetCount), gatingAllHitNeg) + gatingCentral * (microfacetCount-gauss);

#ifdef ALLOW_FLIP
	if (flip)
	{
		float temp = resultPos;
		resultPos = resultNeg;
		resultNeg = temp;
	}
#endif
}

void EvaluateBinomialValueDualGated(float2 randB, float2 randG, float2 successProb, float2 microfacetCount, out float2 resultPos, out float2 resultNeg)
{
	// Fully vectorized implementation

	successProb = saturate(successProb);
#ifdef ALLOW_FLIP
	bool2 flip = successProb > 0.5;
	successProb = lerp(successProb, 1 - successProb, flip);
	// if (flip) successProb = 1 - successProb;
#endif // ALLOW_FLIP

	// Evaluate probabilities for gating
	float2 probAllHitPos = saturate(microfacetCount - 1) * pow(successProb, max(2, microfacetCount));
	float2 probAllHitNeg = saturate(microfacetCount - 1) * OneMinusP_Pow_N(successProb, max(2, microfacetCount));
	float2 probExactOneHitPos = successProb * max(0, 1-abs(1-microfacetCount));
	float2 probExactOneHitNeg = (1-successProb) * max(0, 1-abs(1-microfacetCount));
	float2 probAllOrExactOneHitPos = probAllHitPos + probExactOneHitPos;
	float2 probAllOrExactOneHitNeg = probAllHitNeg + probExactOneHitNeg;

	// Evaluate gatings
	float2 gatingAllHitPos = EvalGating(randB, probAllHitPos);                   // randB < probAllHitPos;
	float2 gatingAllOrOneHitPos = EvalGating(randB, probAllOrExactOneHitPos);    // randB < probAllOrExactOneHitPos
	float2 gatingAllHitNeg = EvalGating(-randB, probAllHitNeg-1);                // randB > (1-probAllHitNeg);
	float2 gatingAllOrOneHitNeg = EvalGating(-randB, probAllOrExactOneHitNeg-1); // randB > (1-probAllOrExactOneHitNeg);
	//float gatingCentral = (microfacetCount > 1) && !(gatingAllOrOneHitPos || gatingAllOrOneHitNeg);
	float2 gatingCentral = (microfacetCount > 1) * (1-(gatingAllOrOneHitPos + gatingAllOrOneHitNeg));

	// Evaluate gaussian approx.
	float2 microfacetCountMean = 1 + (microfacetCount-2) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float2 microfacetCountVar = (microfacetCount-2) * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float2 microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	// Transform to gaussian with mean + stddev.
	float2 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, max(1, microfacetCount-1));

	//float gatingOneHitPos = gatingAllOrOneHitPos - gatingAllHitPos;
	//float gatingOneHitNeg = gatingAllOrOneHitNeg - gatingAllHitNeg;
	// resultPos = gatingOneHitPos * 1 + gatingAllHitPos * max(2, microfacetCount) + gatingCentral * gauss;
	// resultNeg = gatingOneHitNeg * 1 + gatingAllHitNeg * max(2, microfacetCount) + gatingCentral * (microfacetCount-gauss);
	resultPos = gatingAllOrOneHitPos + gatingAllHitPos * max(1, microfacetCount-1) + gatingCentral * gauss;
	resultNeg = gatingAllOrOneHitNeg + gatingAllHitNeg * max(1, microfacetCount-1) + gatingCentral * (max(2, microfacetCount)-gauss);
	//resultPos = lerp(gatingAllOrOneHitPos, max(2, microfacetCount), gatingAllHitPos) + gatingCentral * gauss;
	//resultNeg = lerp(gatingAllOrOneHitNeg, max(2, microfacetCount), gatingAllHitNeg) + gatingCentral * (microfacetCount-gauss);

#ifdef ALLOW_FLIP
	float2 temp = resultPos;
	resultPos = lerp(resultPos, resultNeg, flip);
	resultNeg = lerp(resultNeg, temp, flip);
	/*
	if (flip)
	{
		float2 temp = resultPos;
		resultPos = resultNeg;
		resultNeg = temp;
	}
	*/
#endif
}

void EvaluateBinomialValueDualGated(float4 randB, float4 randG, float4 successProb, float4 microfacetCount, out float4 resultPos, out float4 resultNeg)
{
	// Fully vectorized implementation

	successProb = saturate(successProb);
#ifdef ALLOW_FLIP
	bool4 flip = successProb > 0.5;
	successProb = lerp(successProb, 1 - successProb, flip);
	// if (flip) successProb = 1 - successProb;
#endif // ALLOW_FLIP

	// Evaluate probabilities for gating
	float4 probAllHitPos = saturate(microfacetCount - 1) * pow(successProb, max(2, microfacetCount));
	float4 probAllHitNeg = saturate(microfacetCount - 1) * OneMinusP_Pow_N(successProb, max(2, microfacetCount));
	float4 probExactOneHitPos = successProb * max(0, 1-abs(1-microfacetCount));
	float4 probExactOneHitNeg = (1-successProb) * max(0, 1-abs(1-microfacetCount));
	float4 probAllOrExactOneHitPos = probAllHitPos + probExactOneHitPos;
	float4 probAllOrExactOneHitNeg = probAllHitNeg + probExactOneHitNeg;

	// Evaluate gatings
	float4 gatingAllHitPos = EvalGating(randB, probAllHitPos);                   // randB < probAllHitPos;
	float4 gatingAllOrOneHitPos = EvalGating(randB, probAllOrExactOneHitPos);    // randB < probAllOrExactOneHitPos
	float4 gatingAllHitNeg = EvalGating(-randB, probAllHitNeg-1);                // randB > (1-probAllHitNeg);
	float4 gatingAllOrOneHitNeg = EvalGating(-randB, probAllOrExactOneHitNeg-1); // randB > (1-probAllOrExactOneHitNeg);
	//float gatingCentral = (microfacetCount > 1) && !(gatingAllOrOneHitPos || gatingAllOrOneHitNeg);
	float4 gatingCentral = (microfacetCount > 1) * (1-(gatingAllOrOneHitPos + gatingAllOrOneHitNeg));

	// Evaluate gaussian approx.
	float4 microfacetCountMean = 1 + (microfacetCount-2) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float4 microfacetCountVar = (microfacetCount-2) * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
	float4 microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
	// Transform to gaussian with mean + stddev.
	float4 gauss = randG * microfacetCountStdDev + microfacetCountMean;
	gauss = clamp(gauss, 1, max(1, microfacetCount-1));

	//float gatingOneHitPos = gatingAllOrOneHitPos - gatingAllHitPos;
	//float gatingOneHitNeg = gatingAllOrOneHitNeg - gatingAllHitNeg;
	// resultPos = gatingOneHitPos * 1 + gatingAllHitPos * max(2, microfacetCount) + gatingCentral * gauss;
	// resultNeg = gatingOneHitNeg * 1 + gatingAllHitNeg * max(2, microfacetCount) + gatingCentral * (microfacetCount-gauss);
	resultPos = gatingAllOrOneHitPos + gatingAllHitPos * max(1, microfacetCount-1) + gatingCentral * gauss;
	resultNeg = gatingAllOrOneHitNeg + gatingAllHitNeg * max(1, microfacetCount-1) + gatingCentral * (max(2, microfacetCount)-gauss);
	//resultPos = lerp(gatingAllOrOneHitPos, max(2, microfacetCount), gatingAllHitPos) + gatingCentral * gauss;
	//resultNeg = lerp(gatingAllOrOneHitNeg, max(2, microfacetCount), gatingAllHitNeg) + gatingCentral * (microfacetCount-gauss);

#ifdef ALLOW_FLIP
	float4 temp = resultPos;
	resultPos = lerp(resultPos, resultNeg, flip);
	resultNeg = lerp(resultNeg, temp, flip);
	/*
	if (flip)
	{
		float4 temp = resultPos;
		resultPos = resultNeg;
		resultNeg = temp;
	}
	*/
#endif // ALLOW_FLIP
}

void EvaluateBinomialValueDualGated_Branching(float randB, float randG, float successProb, float microfacetCount, out float resultPos, out float resultNeg)
{
	successProb = saturate(successProb);
	bool flip = successProb > 0.5;
	if (flip) successProb = 1 - successProb;

	//successProb = clamp(successProb, 1e-2, 1-1e-2);
	//microfacetCount = max(0, microfacetCount);
	// TODO merge control flow
	if (microfacetCount <= 1)
	{
		float probAllHitPos = successProb * microfacetCount;
		float probAllHitNeg = (1-successProb) * microfacetCount;
		float probAllHitPosAndNeg = microfacetCount; // <= 1

		float gatingPos = randB < probAllHitPos;
		float gatingNeg = randB < probAllHitNeg;
		float gatingPosAndNeg = randB <= probAllHitPosAndNeg;

		resultPos = gatingPos;
		resultNeg = gatingPosAndNeg - gatingPos;
	}
	else if (microfacetCount <= 2)
	{
		// 5 cases in total!!!
		// TODO optimize!
		float cmfOneHitPos = successProb * (2-microfacetCount);
		float cmfOneHitNeg = (2-microfacetCount);
		float cmfTwoHitPos = cmfOneHitNeg + Sq(successProb) * (microfacetCount-1);
		float cmfTwoHitNeg = cmfTwoHitPos + (Sq(successProb) + Sq(1-successProb)) * (microfacetCount-1);
		float gatingOneHitPos = randB < cmfOneHitPos;
		float gatingOneHitNeg = randB < cmfOneHitNeg;
		float gatingTwoHitPos = randB < cmfTwoHitPos;
		float gatingTwoHitNeg = randB <= cmfTwoHitNeg;
		float gatingOneHitEach = 1.0;
		gatingOneHitEach -= gatingTwoHitNeg;
		gatingTwoHitNeg -= gatingTwoHitPos;
		gatingTwoHitPos -= gatingOneHitNeg;
		gatingOneHitNeg -= gatingOneHitPos;
		resultPos = gatingOneHitEach + gatingOneHitPos + 2*gatingTwoHitPos;
		resultNeg = gatingOneHitEach + gatingOneHitNeg + 2*gatingTwoHitNeg;
	}
	else
	{
		// Compute binomial properties
		float probZeroHitPos = pow(1-successProb, microfacetCount); // Probability of hitting all negative microfacets!
		float probZeroHitNeg = pow(successProb, microfacetCount);   // Probability of hitting all positive microfacets!
		//float probOneHitEach = 1 - probZeroHitPos - probZeroHitNeg; // Probability of hitting at leat one negative microfacet (given that we have hit at least one positive microfacet).
		float cmfZeroHitPos = probZeroHitPos;
		float cmfOneHitEach = 1 - probZeroHitNeg; // probZeroHitPos + probOneHitEach;

		float microfacetCountMean = 1 + (microfacetCount-2) * successProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
		float microfacetCountVar = (microfacetCount-2) * successProb * (1-successProb); // Standard deviation of number of hits in the footprint given already one hit
		float microfacetCountStdDev = sqrt(max(microfacetCountVar, 0));
		float binomialSmoothWidthPos = 0; //0.1 * clamp(probOneHitPos * 10, 0.0, 1.0) * clamp((1.0 - probOneHitPos) * 10, 0.0, 1.0);
		float binomialSmoothWidthNeg = 0; //0.1 * clamp(probOneHitNeg * 10, 0.0, 1.0) * clamp((1.0 - probOneHitNeg) * 10, 0.0, 1.0);

		// 3 cases
		// TODO optimize!
		float gatingZeroHitPos = randB < cmfZeroHitPos;
		float gatingOneHitEach = randB <= cmfOneHitEach;
		float gatingZeroHitNeg = 1.0;
		gatingZeroHitNeg -= gatingOneHitEach;
		gatingOneHitEach -= gatingZeroHitPos;

		// Transform to gaussian with mean + stddev.
		float gauss = randG * microfacetCountStdDev + microfacetCountMean;
		// NOTE: We have adjusted the mean such that we don't need to adjust the clamp limits!
		// YES WE DO NEED TO ADJUST CLAMP LIMITS! at least one microfacet is left and at least one microfacet is right!
		gauss = clamp(gauss, 1, max(1, microfacetCount-1));
		// Combine the three gating cases...
		resultPos = /* gatingZeroHitPos * 0 + */ gatingOneHitEach * gauss + gatingZeroHitNeg * microfacetCount;
		resultNeg = microfacetCount-resultPos;
		//resultNeg = /* gatingZeroHitNeg * 0 + */ gatingOneHitEach * (microfacetCount - gauss) + gatingZeroHitPos * microfacetCount;
	}

	if (flip)
	{
		float temp = resultPos;
		resultPos = resultNeg;
		resultNeg = temp;
	}
}



#endif // GLINT_BINOMIAL_HLSL
