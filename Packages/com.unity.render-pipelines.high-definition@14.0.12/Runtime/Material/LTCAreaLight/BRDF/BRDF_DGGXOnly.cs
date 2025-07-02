using System;
using UnityEngine.Rendering;

namespace UnityEngine.Rendering.HighDefinition.LTC
{
    /// <summary>
    /// GGX implementation of the BRDF interface
    /// </summary>
    class BRDF_DGGXOnly : IBRDF
    {
        public double Eval(ref Vector3 _tsView, ref Vector3 _tsLight, float _alpha, out double _pdf)
        {
            if (_tsView.z <= 0)
            {
                _pdf = 0;
                return 0;
            }

            // D
            Vector3 H = _tsView + _tsLight;
            float lengthH = H.magnitude;
            if (lengthH > 1e-8f)
                H = H / lengthH;
            else
                H = new Vector3(0, 0, 1);

            double slopex = H.x / H.z;
            double slopey = H.y / H.z;
            double D = 1.0 / (1.0 + (slopex * slopex + slopey * slopey) / _alpha / _alpha);
            D = D * D;
            D = D / (Math.PI * _alpha * _alpha * H.z * H.z * H.z * H.z);

            double HdotV = Vector3.Dot(_tsView, H);

            // We want to simply integrate \int_{\domainH} D(\dirH) \dx{\dirH} over the light source
            // Since the domain is \domainL here, we perform a change of variables and we get
            // \int_{\domainL} D(\dirH) / (4*HdotV) \dx{\dirL}
            double res = D / (4.0 * HdotV);

            // pdf = D(H) * (N.H) / (4 * (L.H))
            _pdf = Math.Abs(D * H.z / (4.0 * HdotV));

            return res;
        }

        public void GetSamplingDirection(ref Vector3 _tsView, float _alpha, float _U1, float _U2, ref Vector3 _direction)
        {
            float phi = 2.0f * Mathf.PI * _U1;
            float r = _alpha * Mathf.Sqrt(_U2 / (1.0f - _U2));
            Vector3 H = new Vector3(r * Mathf.Cos(phi), r * Mathf.Sin(phi), 1.0f).normalized;
            _direction = -_tsView + 2.0f * H * Vector3.Dot(H, _tsView);
        }

        public LTCLightingModel GetLightingModel()
        {
            return LTCLightingModel.DGGXOnly;
        }
    }
}
