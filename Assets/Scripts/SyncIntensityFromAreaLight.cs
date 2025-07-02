using UnityEngine;
using UnityEngine.Rendering.HighDefinition;

internal static class MaterialExtension
{
    public static void UpdateEmissiveColorFromIntensityAndEmissiveColorLDR(this Material material)
    {
        const string kEmissiveColorLDR = "_EmissiveColorLDR";
        const string kEmissiveColor = "_EmissiveColor";
        const string kEmissiveIntensity = "_EmissiveIntensity";

        if (material.HasProperty(kEmissiveColorLDR) && material.HasProperty(kEmissiveIntensity) && material.HasProperty(kEmissiveColor))
        {
            // Important: The color picker for kEmissiveColorLDR is LDR and in sRGB color space but Unity don't perform any color space conversion in the color
            // picker BUT only when sending the color data to the shader... So as we are doing our own calculation here in C#, we must do the conversion ourselves.
            Color emissiveColorLDR = material.GetColor(kEmissiveColorLDR);
            Color emissiveColorLDRLinear = new Color(Mathf.GammaToLinearSpace(emissiveColorLDR.r), Mathf.GammaToLinearSpace(emissiveColorLDR.g), Mathf.GammaToLinearSpace(emissiveColorLDR.b));
            material.SetColor(kEmissiveColor, emissiveColorLDRLinear * material.GetFloat(kEmissiveIntensity));
        }
    }
}


[ExecuteInEditMode]
public class SyncIntensityFromAreaLight : MonoBehaviour
{
    void Update()
    {
        var light = GetComponent<Light>();
        var intensity = light.intensity; // In Nits

        var renderer = GetComponentInChildren<MeshRenderer>();
        var material = renderer.sharedMaterial;

        material.SetInt("_EmissiveIntensityUnit", 0);
        material.SetInt("_UseEmissiveIntensity", 1);
        material.SetFloat("_EmissiveIntensity", intensity);
        material.UpdateEmissiveColorFromIntensityAndEmissiveColorLDR();
    }
}
