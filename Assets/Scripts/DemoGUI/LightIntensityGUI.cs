using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class LightIntensityGUI : MonoBehaviour
{
    public HDAdditionalLightData target;

    public Slider slider;

    public TextMeshProUGUI label;
    public string format;

    public void Start()
    {
        if (slider != null)
        {
            if (target != null)
            {
                slider.value = Mathf.Log10(target.intensity);
            }

            slider.onValueChanged.AddListener(delegate { SetValue(slider.value); });
        }
        if (label != null)
        {
            label.text = string.Format(format, Mathf.Log10(target.intensity));
        }
    }

    public void SetValue(float value)
    {
        if (target != null)
        {
            target.intensity = Mathf.Pow(10.0f, value);
        }
        if (label != null)
        {
            label.text = string.Format(format, value);
        }
    }
}
