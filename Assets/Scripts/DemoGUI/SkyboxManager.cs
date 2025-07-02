using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System.Collections.Generic;
using TMPro;

public class SkyboxManager : MonoBehaviour
{
    public Volume volume;
    public VolumeProfile[] volumeProfiles;
    public string[] volumeProfileNames;

    public int currentIndex = 0;

    public TMP_Dropdown dropdownMenu;

    public void Start()
    {
        if (dropdownMenu != null)
        {
            dropdownMenu.options.Clear();
            for (int i = 0; i < volumeProfileNames.Length; ++i)
            {
                dropdownMenu.options.Add(new TMP_Dropdown.OptionData(volumeProfileNames[i]));
            }

            dropdownMenu.onValueChanged.AddListener(delegate { SetCurrentIndex(dropdownMenu.value); });
        }
        AssignVolume();
    }

    public void SetCurrentIndex(int value)
    {
        currentIndex = value;
        AssignVolume();
    }

    private void AssignVolume()
    {
        if (currentIndex >= volumeProfiles.Length)
        {
            Debug.LogError($"Current sky index out of bounds {currentIndex} length={volumeProfiles.Length}");
            return;
        }
        if (volume == null)
        {
            Debug.LogError("Volume not assigned!");
            return;
        }
        volume.sharedProfile = volumeProfiles[currentIndex];
    }
}
