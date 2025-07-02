using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System.Collections.Generic;
using TMPro;

public class MaterialManager : MonoBehaviour
{
    public Renderer[] renderers;
    public Material[] materials;
    public string[] materialNames;

    public int currentIndex = 0;

    public TMP_Dropdown dropdownMenu;

    public void Start()
    {
        if (dropdownMenu != null)
        {
            dropdownMenu.options.Clear();
            for (int i = 0; i < materialNames.Length; ++i)
            {
                dropdownMenu.options.Add(new TMP_Dropdown.OptionData(materialNames[i]));
            }

            dropdownMenu.onValueChanged.AddListener(delegate { SetCurrentIndex(dropdownMenu.value); });
        }

        AssignMaterial();
    }

    public void SetCurrentIndex(int value)
    {
        currentIndex = value;
        AssignMaterial();
    }

    public Material GetMaterial()
    {
        if (currentIndex >= materials.Length)
        {
            Debug.LogError($"Current material index out of bounds {currentIndex} length={materials.Length}");
            return null;
        }
        return materials[currentIndex];
    }

    private void AssignMaterial()
    {
        if (currentIndex >= materials.Length)
        {
            Debug.LogError($"Current material index out of bounds {currentIndex} length={materials.Length}");
            return;
        }
        foreach (var renderer in this.renderers)
        {
            renderer.sharedMaterial = materials[currentIndex];
        }
    }
}
