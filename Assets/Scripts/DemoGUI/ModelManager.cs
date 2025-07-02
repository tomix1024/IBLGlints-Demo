using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.HighDefinition;
using System.Collections.Generic;
using TMPro;

public class ModelManager : MonoBehaviour
{
    public GameObject[] gameObjects;
    public string[] gameObjectNames;

    public int currentIndex = 0;

    public TMP_Dropdown dropdownMenu;

    public void Start()
    {
        if (dropdownMenu != null)
        {
            dropdownMenu.options.Clear();
            for (int i = 0; i < gameObjectNames.Length; ++i)
            {
                dropdownMenu.options.Add(new TMP_Dropdown.OptionData(gameObjectNames[i]));
            }

            dropdownMenu.onValueChanged.AddListener(delegate { SetCurrentIndex(dropdownMenu.value); });
        }

        UpdateActive();
    }

    public void SetCurrentIndex(int value)
    {
        currentIndex = value;
        UpdateActive();
    }

    private void UpdateActive()
    {
        if (currentIndex >= gameObjects.Length)
        {
            Debug.LogError($"Current object index out of bounds {currentIndex} length={gameObjects.Length}");
            return;
        }
        for (int i = 0; i < gameObjects.Length; ++i)
        {
            gameObjects[i].SetActive(currentIndex == i);
        }
    }
}
