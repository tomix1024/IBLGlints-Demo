Real-time Image-based Lighting of Glints (Demo)
===============================================


Overview
--------

We have introduced an extension to real-time image-based lighting at EGSR 2025.
This repository contains a Unity3D project that demonstrates real-time rendering of glints under image-based lighting [KK25], area lights [KK24] and directional lights [DB23] (modified based on area-light theory).
We provide a modified version of Unity3D's High Definition Render Pipeline containing our implementation [including parts of previous work](https://thomasdeliot.wixsite.com/blog/single-post/hpg23-real-time-rendering-of-glinty-appearance-using-distributed-binomial-laws-on-anisotropic-grids).

> - [KK25] Tom Kneiphof and Reinhard Klein. "Real-time Image-based Lighting of Glints." Computer Graphics Forum. Vol. 44. No. 4. 2025.
> - [KK24] Tom Kneiphof and Reinhard Klein. "Real-Time Rendering of Glints in the Presence of Area Lights." Pacific Graphics 2024.
> - [DB23] Thomas Deliot and Laurent Belcour. "Real‐Time Rendering of Glinty Appearances using Distributed Binomial Laws on Anisotropic Grids." Computer Graphics Forum. Vol. 42. No. 8. 2023.


Project Structure
-----------------

A copy of the HDRP and core render pipelines are found in the `Packages/` directory.
Our implementation modifies the `HDRP/Lit` shader in the HDRP.
The relevant entry points for the glint evaluation are found in `Packages/com.unity.render-pipelines.high-definition@14.0.10/Runtime/Material/Lit/Lit.hlsl` in the functions
- `IndirectLighting EvaluateBSDF_Env(...)`
- `DirectLighting EvaluateBSDF_Rect(...)`
- `CBSDF EvaluateBSDF(...)`

The binomial/multinomial sampling is then implemented in the files `Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Glints*.hlsl`.

The Unity Project itself contains the demo scene that can be built into an executable.


Precompiled Binaries
--------------------

Precompiled binaries are available for [Linux](https://cg.cs.uni-bonn.de/backend/v1/files/code/IBLGlints-Demo/IBLGlints-Demo-Linux.tar.gz) and [Windows](https://cg.cs.uni-bonn.de/backend/v1/files/code/IBLGlints-Demo/IBLGlints-Demo-Windows.zip).

| Filename                    | sha256sum                                                        |
|-----------------------------|------------------------------------------------------------------|
| [IBLGlints-Demo-Linux.tar.gz](https://cg.cs.uni-bonn.de/backend/v1/files/code/IBLGlints-Demo/IBLGlints-Demo-Linux.tar.gz) | e87a87201b7fe76291fa95419fa3db6e924f013e676f3ec09c79b695c5041dd3 |
| [IBLGlints-Demo-Windows.zip](https://cg.cs.uni-bonn.de/backend/v1/files/code/IBLGlints-Demo/IBLGlints-Demo-Windows.zip)  | 052bc4170e70606fd19a579c737a8f4e2250a1ab6e4fece951b5a8c7df8d40c7 |

Citation
--------

If you are using this code in academic research, please cite our paper.
The BibTeX entry is
```bibtex
@article{kneiphof2025real,
  title={Real-time Image-based Lighting of Glints},
  author={Kneiphof, Tom and Klein, Reinhard},
  journal={Computer Graphics Forum},
  year = {2025},
  publisher = {The Eurographics Association and John Wiley & Sons Ltd.},
  ISSN = {1467-8659},
  DOI = {10.1111/cgf.70175}
}
```
