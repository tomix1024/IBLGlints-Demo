using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class SkyRotationGUI : MonoBehaviour
{
    public Volume target;

    public Slider slider;

    public TextMeshProUGUI label;
    public string format;

    private VolumeProfile profile;
    private SkySettings skySettings;
    private float rotationValue;

    private void UpdateSkySettings()
    {
        if (target.sharedProfile != profile)
        {
            profile = target.sharedProfile;
            profile.TryGetSubclassOf(typeof(SkySettings), out skySettings);
            rotationValue = skySettings.rotation.value;
            slider.value = rotationValue;
            label.text = string.Format(format, rotationValue);
        }
    }

    public void Start()
    {
        UpdateSkySettings();

        if (slider != null)
        {
            slider.value = rotationValue;
            slider.onValueChanged.AddListener(delegate { SetRotationValue(slider.value); });
        }
    }

    public void Update()
    {
        UpdateSkySettings();
    }

    public void SetRotationValue(float value)
    {
        rotationValue = value;

        skySettings.rotation.overrideState = true;
        skySettings.rotation.value = rotationValue;

        label.text = string.Format(format, rotationValue);
    }
}
