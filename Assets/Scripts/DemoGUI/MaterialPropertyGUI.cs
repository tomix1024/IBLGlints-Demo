using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class MaterialPropertyGUI : MonoBehaviour
{
    public Material[] materials;

    public string propertyName;

    public Slider slider;

    public TextMeshProUGUI label;
    public string format;

    public bool logarithmic = false;
    public bool decodeLogSinSunAngle = false;

    public void Start()
    {
        if (slider != null)
        {
            Material material = materials[0];
            if (material != null)
            {
                float value = material.GetFloat(propertyName);
                if (logarithmic)
                {
                    Debug.Log($"before log {value}");
                    value = Mathf.Log(value, 2.0f);
                    Debug.Log($"after log {value}");
                }
                slider.value = value;
            }

            slider.onValueChanged.AddListener(delegate { SetValue(slider.value); });
        }

        if (label != null)
        {
            Material material = materials[0];
            if (material != null)
            {
                float value = material.GetFloat(propertyName);
                if (decodeLogSinSunAngle)
                    value = Mathf.Rad2Deg * Mathf.Asin(Mathf.Exp(value));
                string text = string.Format(format, value); //value.ToString("0.00");
                label.text = text;
            }
        }
    }

    public void SetValue(float value)
    {
        if (logarithmic)
        {
            Debug.Log($"before pow {value}");
            value = Mathf.Pow(2.0f, value);
            Debug.Log($"after pow {value}");
        }
        foreach (var material in materials)
        {
            if (material.HasProperty(propertyName))
            {
                material.SetFloat(propertyName, value);
                if (decodeLogSinSunAngle)
                    ValidateMaterial(material);
            }
        }

        if (label != null)
        {
            if (decodeLogSinSunAngle)
                value = Mathf.Rad2Deg * Mathf.Asin(Mathf.Exp(value));
            string text = string.Format(format, value); //value.ToString("0.00");
            label.text = text;
        }
    }

    const string kLogSinSunAngle = "_LogSinSunAngle";
    const string kSunSolidAngle = "_SunSolidAngle";
    private static void ValidateMaterial(Material material)
    {
        if (material.HasProperty(kLogSinSunAngle) && material.HasProperty(kSunSolidAngle))
        {
            float logSinSunAngle = material.GetFloat(kLogSinSunAngle);
            float sinGammaSq = Mathf.Exp(2*logSinSunAngle);
            float cosGamma = Mathf.Sqrt(1 - sinGammaSq);
            float Al = 2*Mathf.PI * (1 - cosGamma);
            material.SetFloat(kSunSolidAngle, Al);
        }
    }
}
