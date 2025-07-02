using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class ExposureGUI : MonoBehaviour
{
    public Volume target;

    public Slider slider;

    public TextMeshProUGUI label;
    public string format;

    private float exposureValue;

    public void Start()
    {
        if (target != null)
        {
            VolumeProfile profile = target.sharedProfile;
            Exposure exposureSettings;
            if (profile.TryGet<Exposure>(out exposureSettings))
            {
                exposureValue = -exposureSettings.fixedExposure.value;
            }
        }
        if (slider != null)
        {

            slider.value = exposureValue;
            slider.onValueChanged.AddListener(delegate { SetValue(slider.value); });
        }
        if (label != null)
        {
            label.text = string.Format(format, exposureValue);
        }
    }

    public void SetValue(float value)
    {
        exposureValue = value;
        if (target != null)
        {
            VolumeProfile profile = target.sharedProfile;
            Exposure exposureSettings;
            if (profile.TryGet<Exposure>(out exposureSettings))
            {
                exposureSettings.fixedExposure.value = -exposureValue;
            }
        }
        if (label != null)
        {
            label.text = string.Format(format, exposureValue);
        }
    }
}
