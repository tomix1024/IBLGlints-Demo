using UnityEngine;
using UnityEngine.UI;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System;
using System.Collections.Generic;
using TMPro;

public class MaterialEnumGUI : MonoBehaviour
{
    public Material[] materials;

    public string[] keys;
    public int[] values;

    public string propertyName;

    public TMP_Dropdown dropdownMenu;
    public int currentIndex = 0;

    public void Start()
    {
        if (dropdownMenu != null)
        {
            dropdownMenu.options.Clear();
            for (int i = 0; i < keys.Length; ++i)
            {
                dropdownMenu.options.Add(new TMP_Dropdown.OptionData(keys[i]));
            }

            dropdownMenu.onValueChanged.AddListener(delegate { SetCurrentIndex(dropdownMenu.value); });
        }
    }

    public void SetCurrentIndex(int value)
    {
        currentIndex = value;
        AssignMaterial();
    }

    private void AssignMaterial()
    {
        if (currentIndex >= values.Length)
        {
            Debug.LogError($"Current index out of bounds {currentIndex} length={values.Length}");
            return;
        }

        foreach (var material in materials)
        {
            material.SetInt(propertyName, values[currentIndex]);
            ValidateMaterial(material);
        }
    }

    const string kGlintIBLVariant = "_GlintIBLVariant";
    private static void ValidateMaterial(Material material)
    {
        if (material.HasProperty(kGlintIBLVariant))
        {
            int iblVariant = material.GetInt(kGlintIBLVariant);
            CoreUtils.SetKeyword(material, "_GLINTS_IBL_VARIANT_DISCRETE4", iblVariant == 0);
            CoreUtils.SetKeyword(material, "_GLINTS_IBL_VARIANT_DISCRETE8", iblVariant == 1);
            CoreUtils.SetKeyword(material, "_GLINTS_IBL_VARIANT_WANG2020", iblVariant == 2);
        }
    }
}
