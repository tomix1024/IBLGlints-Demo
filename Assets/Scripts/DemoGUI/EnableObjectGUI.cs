using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class EnableObjectGUI : MonoBehaviour
{
    public GameObject target;

    public Toggle checkbox;

    public void Start()
    {
        if (checkbox != null)
        {
            if (target != null)
            {
                checkbox.isOn = target.activeSelf;
            }

            checkbox.onValueChanged.AddListener(delegate { SetValue(checkbox.isOn); });
        }
    }

    public void SetValue(bool value)
    {
        if (checkbox == null)
            return;
        if (target == null)
            return;

        target.SetActive(value);
    }
}
