using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class AreaLightGUI : MonoBehaviour
{
    public HDAdditionalLightData target;

    public Slider sizeSlider;
    public Slider ratioSlider;

    public TextMeshProUGUI sizeLabel;
    public TextMeshProUGUI ratioLabel;
    public string format;

    private float logSize;
    private float logRatio;

    public void Start()
    {
        if (target != null)
        {
            logSize = Mathf.Log(target.shapeWidth * target.shapeHeight, 2.0f);
            logRatio = Mathf.Log(target.shapeWidth / target.shapeHeight, 2.0f);
        }

        if (sizeSlider != null)
        {
            sizeSlider.value = logSize;
            sizeSlider.onValueChanged.AddListener(delegate { SetLogSize(sizeSlider.value); });
        }
        if (sizeLabel != null)
        {
            sizeLabel.text = string.Format(format, logSize);
        }
        if (ratioSlider != null)
        {
            ratioSlider.value = logRatio;
            ratioSlider.onValueChanged.AddListener(delegate { SetLogRatio(ratioSlider.value); });
        }
        if (ratioLabel != null)
        {
            ratioLabel.text = string.Format(format, logRatio);
        }
    }

    public void SetLogSize(float value)
    {
        logSize = value;
        if (target != null)
        {
            target.shapeWidth = Mathf.Pow(2.0f, (logSize+logRatio)/2);
            target.shapeHeight = Mathf.Pow(2.0f, (logSize-logRatio)/2);
        }
        if (sizeLabel != null)
        {
            sizeLabel.text = string.Format(format, logSize);
        }
    }
    public void SetLogRatio(float value)
    {
        logRatio = value;
        if (target != null)
        {
            target.shapeWidth = Mathf.Pow(2.0f, (logSize+logRatio)/2);
            target.shapeHeight = Mathf.Pow(2.0f, (logSize-logRatio)/2);
        }
        if (ratioLabel != null)
        {
            ratioLabel.text = string.Format(format, logRatio);
        }
    }
}
