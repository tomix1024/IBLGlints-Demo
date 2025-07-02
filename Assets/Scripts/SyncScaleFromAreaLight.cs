using UnityEngine;
using UnityEngine.Rendering.HighDefinition;

[ExecuteInEditMode]
public class SyncScaleFromAreaLight : MonoBehaviour
{
    void Update()
    {
        var light = GetComponent<HDAdditionalLightData>();
        var width = light.shapeWidth;
        var height = light.shapeHeight;
        if (transform.localScale.x != width || transform.localScale.y != height)
            transform.localScale = new Vector3(width, height, 1);
    }
}
