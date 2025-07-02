using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class TransformGUI : MonoBehaviour
{
    public Transform target;

    public Slider polarSlider;
    public Slider azimuthSlider;

    public TextMeshProUGUI polarLabel;
    public TextMeshProUGUI azimuthLabel;
    public string format;

    public void Start()
    {
        // x = polar (deg)
        // y = azimuth (deg)
        if (polarSlider != null)
        {
            if (target != null)
            {
                polarSlider.value = target.localEulerAngles.z;
            }

            polarSlider.onValueChanged.AddListener(delegate { SetPolarValue(polarSlider.value); });
        }
        if (azimuthSlider != null)
        {
            if (target != null)
            {
                azimuthSlider.value = target.localEulerAngles.y;
            }

            azimuthSlider.onValueChanged.AddListener(delegate { SetAzimuthValue(azimuthSlider.value); });
        }
        if (polarLabel != null)
        {
            polarLabel.text = string.Format(format, target.localEulerAngles.z);
        }
        if (azimuthLabel != null)
        {
            azimuthLabel.text = string.Format(format, target.localEulerAngles.y);
        }
    }

    public void SetPolarValue(float value)
    {
        if (target != null)
        {
            target.localEulerAngles = new Vector3(target.localEulerAngles.x, target.localEulerAngles.y, value);
        }
        if (polarLabel != null)
        {
            polarLabel.text = string.Format(format, value);
        }
    }
    public void SetAzimuthValue(float value)
    {
        if (target != null)
        {
            target.localEulerAngles = new Vector3(target.localEulerAngles.x, value, target.localEulerAngles.z);
        }
        if (azimuthLabel != null)
        {
            azimuthLabel.text = string.Format(format, value);
        }
    }
}
