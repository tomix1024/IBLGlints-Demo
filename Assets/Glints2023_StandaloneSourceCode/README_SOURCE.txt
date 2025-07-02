The source code contains two parts:

- An HLSL file to include in your rendering shader, that contains the Glinty NDF function to call during your light loop.
- A C# script combined with a shader file to generate the small precomputed texture of packed random numbers that we use to speed up performance. Since we compute this texture on the GPU, it takes only a few milliseconds and only needs to be run once. We provide this as a Unity specific implementation, but the code is short and easy to port to other frameworks.

To use the Glinty NDF, simply replace the NDF evaluation in your code with a call to the SampleGlints2023NDF method. The input parameters are the following:
- localHalfVector: the half vector in tangent space
- targetNDF: the target NDF energy D(h)
- maxNDF: the maximym possible NDF energy D(0) (when the local half vector's slope is null)
- uv: the UV texture coordinates
- duvdx, duvdy: the screen-space derivatives of UV

The include file additionally defines the following constant parameters that need to be sent to the shader:
_Glint2023NoiseMap: the precomputed texture of random numbers to bind
_Glint2023NoiseMapSize: the width of the precomputed texture of random numbers (the texture is assumed to be square)
_ScreenSpaceScale: the minimal allowed size of a glint in pixels. Due to Shannon's sampling theorem, using 1 here results in aliasing. We recommend using 1.5
_LogMicrofacetDensity: the (log) amount of microfacets per unit surface texture space area 
_MicrofacetRoughness: the ad-hoc R ratio parameter in our paper, simulating a microfacet roughness parameter
_DensityRandomization: the amount of randomization of glint densities across the surface