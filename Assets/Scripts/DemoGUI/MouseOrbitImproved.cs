using UnityEngine;
using System.Collections;
using UnityEngine.EventSystems;

[AddComponentMenu("Camera-Control/Mouse Orbit with zoom")]
public class MouseOrbitImproved : MonoBehaviour {
 
    public Transform target;
    public float logDistance = 0.5f;
    public float logDistanceTarget = 0.5f;

    public float xSpeed = 5.0f;
    public float ySpeed = 5.0f;
    public float distanceSpeed = 5.0f;
 
    public float yMinLimit = -20f;
    public float yMaxLimit = 80f;
 
    public float logDistanceMin = -1.0f;
    public float logDistanceMax = 5.0f;
    public float distanceInterpolationMass = 100.0f;

    private Rigidbody rigidbody;
 
    private bool onGUI = false;

    float x = 0.0f;
    float y = 0.0f;
 
    // Use this for initialization
    void Start () 
    {
        Vector3 angles = transform.eulerAngles;
        x = angles.y;
        y = angles.x;
 
        rigidbody = GetComponent<Rigidbody>();
 
        // Make the rigid body not change rotation
        if (rigidbody != null)
        {
            rigidbody.freezeRotation = true;
        }
    }
 
    void FixedUpdate()
    {
        // C.f. exponential moving average, optical depth integration, etc.
        float alpha = Mathf.Exp(-Time.fixedDeltaTime * distanceInterpolationMass);
        logDistance = Mathf.Lerp(logDistance, logDistanceTarget, alpha);
    }

    void LateUpdate()
    {
        if (target) 
        {
            if (Input.GetMouseButtonDown(0))
                onGUI = EventSystem.current.IsPointerOverGameObject();

            if (Input.GetMouseButton(0) && !onGUI)
            {
                x += Input.GetAxis("Mouse X") * xSpeed;
                y -= Input.GetAxis("Mouse Y") * ySpeed;
            }
 
            y = ClampAngle(y, yMinLimit, yMaxLimit);
 
            Quaternion rotation = Quaternion.Euler(y, x, 0);

            if (!onGUI)
            {
                logDistanceTarget = Mathf.Clamp(logDistanceTarget - Input.GetAxis("Mouse ScrollWheel")*distanceSpeed, logDistanceMin, logDistanceMax);
            }

            // FixedUpdate logDistanceTarget -> logDistance....

            Vector3 negDistance = new Vector3(0.0f, 0.0f, -Mathf.Exp(logDistance));
            Vector3 position = rotation * negDistance + target.position;
 
            transform.rotation = rotation;
            transform.position = position;
        }
    }
 
    public static float ClampAngle(float angle, float min, float max)
    {
        if (angle < -360F)
            angle += 360F;
        if (angle > 360F)
            angle -= 360F;
        return Mathf.Clamp(angle, min, max);
    }
}
