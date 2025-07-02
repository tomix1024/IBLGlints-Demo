using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[ExecuteInEditMode]
public class Glints2023: MonoBehaviour
{
	[Header("Interface")]
	public bool resetButton = false;

	[Header("Noise Settings")]
	public int noiseTexSize = 512;
	public Texture2D glintNoiseTex;
	private Material glintNoiseInitMaterial;


	void Start()
	{
		ResetEverything();
	}

	void OnEnable()
	{
		ResetEverything();
	}

	void OnDisable()
	{

	}

	void Update()
	{
		if(resetButton == true)
		{
			ResetEverything();
		}

#if UNITY_EDITOR
		Shader.SetGlobalTexture("_Glint2023NoiseMap", glintNoiseTex);
		Shader.SetGlobalInt("_Glint2023NoiseMapSize", noiseTexSize);
#endif
	}



	private void ResetEverything()
	{
		resetButton = false;
		OnDisable();

#if UNITY_EDITOR
		glintNoiseInitMaterial = new Material(Shader.Find("Custom/Glints2023NoiseInit"));
		GenerateGlintNoiseTex();
		glintNoiseTex = AssetDatabase.LoadAssetAtPath<Texture2D>("Assets/glint2023Noise.asset");
#endif
		Shader.SetGlobalTexture("_Glint2023NoiseMap", glintNoiseTex);
		Shader.SetGlobalInt("_Glint2023NoiseMapSize", noiseTexSize);
	}

	float InvCDF(float U, float mu, float sigma)
	{
		float x = sigma * Mathf.Sqrt(2.0f) * ErfInv(2.0f * U - 1.0f) + mu;
		return x;
	}

	float ErfInv(float x)
	{
		float w, p;
		w = -Mathf.Log((1.0f - x) * (1.0f + x));
		if (w < 5.000000f)
		{
			w = w - 2.500000f;
			p = 2.81022636e-08f;
			p = 3.43273939e-07f + p * w;
			p = -3.5233877e-06f + p * w;
			p = -4.39150654e-06f + p * w;
			p = 0.00021858087f + p * w;
			p = -0.00125372503f + p * w;
			p = -0.00417768164f + p * w;
			p = 0.246640727f + p * w;
			p = 1.50140941f + p * w;
		}
		else
		{
			w = Mathf.Sqrt(w) - 3.000000f;
			p = -0.000200214257f;
			p = 0.000100950558f + p * w;
			p = 0.00134934322f + p * w;
			p = -0.00367342844f + p * w;
			p = 0.00573950773f + p * w;
			p = -0.0076224613f + p * w;
			p = 0.00943887047f + p * w;
			p = 1.00167406f + p * w;
			p = 2.83297682f + p * w;
		}
		return p * x;
	}

	private void GenerateGlintNoiseTex()
	{
		// Generate noise
		RenderTexture renderTex = new RenderTexture(noiseTexSize, noiseTexSize, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Linear);
		renderTex.enableRandomWrite = true;
		renderTex.useMipMap = false;
		renderTex.autoGenerateMips = false;
		renderTex.Create();
		glintNoiseInitMaterial.SetVector("_FrameSize", new Vector4(noiseTexSize, noiseTexSize, 0, 0));
		glintNoiseInitMaterial.SetInt("_Seed", (int)(Random.value * 100));
		Graphics.Blit(null, renderTex, glintNoiseInitMaterial);

		// Apply to texture
		glintNoiseTex = new Texture2D(noiseTexSize, noiseTexSize, TextureFormat.RGBAFloat, false, true);
		glintNoiseTex.name = "NoiseMap";
		glintNoiseTex.filterMode = FilterMode.Point;
		glintNoiseTex.anisoLevel = 1;
		glintNoiseTex.wrapMode = TextureWrapMode.Repeat;
		glintNoiseTex.ReadPixels(new Rect(0, 0, noiseTexSize, noiseTexSize), 0, 0);
		glintNoiseTex.Apply(false);

#if UNITY_EDITOR
		AssetDatabase.CreateAsset(glintNoiseTex, "Assets/glint2023Noise.asset");
		AssetDatabase.SaveAssets();
		AssetDatabase.Refresh();
#endif

		RenderTexture.active = null;
		renderTex.Release();
	}
}
