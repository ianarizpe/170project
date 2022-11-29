// simple_scene.fs

#version 410 core

uniform vec2 iRes;
uniform float iTime;
out vec4 fragColor;

#define MAX_STEPS 1000
#define MAX_DIST 1000
#define SURF_DIST 0.01

#define M_PI 3.1415926535897932384626433832795

// SDF for torus 
float sdfTorus(vec3 p, vec2 r)
{
  vec2 q = vec2(p.y, length(p.xz) - r.x);
  float d = length(q) - r.y;
  return d;
}

float sdPlane( vec3 p, vec3 n, float h )
{
  // n must be normalized
  return dot(p,n) + h ;
}

float sdCone( in vec3 p, in vec2 c, float h )
{
  // c is the sin/cos of the angle, h is height
  // Alternatively pass q instead of (c,h),
  // which is the point at the base in 2D
  vec2 q = h*vec2(c.x/c.y,-1.0);
    
  vec2 w = vec2( length(p.xz), p.y );
  vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
  float k = sign( q.y );
  float d = min(dot( a, a ),dot(b, b));
  float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
  return sqrt(d)*sign(s);
}

float getDist(vec3 p)
{
  // distance to plane
  float pd = sdPlane(p, vec3(0, 0, 1), 0);

  vec2 r = vec2(1, 0.25);
  float td = sdfTorus(p, r);

  // (sin(t), cos(t))
  vec2 c = vec2(0.49999999999999994, 0.8660254037844387);
  float cd = sdCone(p, c, 2);

  // return min 
  return cd; //min(pd, td);
}

float rayMarch(vec3 ro, vec3 rd)
{
  float d0 = 0.0;

  for (int i = 0; i < MAX_STEPS; i++) {
    vec3 p = ro + d0*rd;
    float ds = getDist(p);
    d0 += ds;
    if (d0 > MAX_DIST || ds < SURF_DIST)
      break;
  }

  return d0;    
}

/* 
  getNormal()

  Computes surface normal using the gradient formula.
*/
vec3 getNormal(vec3 p)
{
  // get distance 
  float d = getDist(p);
  
  // define epsilion
  vec2 e = vec2(0.01, 0.0);

  vec3 n = vec3 (getDist(p + e.xyy) - getDist(p - e.xyy),
                 getDist(p + e.yxy) - getDist(p - e.yxy),
                 getDist(p + e.yyx) - getDist(p - e.yyx));
  // normalize vector 
  return normalize(n);
}

float getLight(vec3 p)
{
  // light position
  vec3 lightPos = vec3(5, 5, 0);

  // animate light pos
  vec2 lr = 5.0 * vec2(sin(iTime), cos(iTime));
  lightPos.xz += lr;

  vec3 l = normalize(lightPos - p);
  vec3 n = getNormal(p);
  // clamp to avoid negative values 
  float dif = clamp(dot(l, n), 0.0, 1.0);

  /*
  // add shadow 
  // start marching a bit high above the surface 
  float d = rayMarch(p + n*SURF_DIST*2, l);
  if (d < length(lightPos - p)) {
    dif *= 0.1;
  }*/

  return dif;
}


mat4 LookAt(vec3 eye, vec3 at, vec3 up)
{
  vec3 Z = normalize(at - eye);    
  vec3 X = normalize(cross(Z, up));
  vec3 Y = cross(X, Z);

  Z.xyz = -Z.xyz;

  mat4 T = mat4(
    vec4(1.0, 0.0, 0.0, 0.0),
    vec4(0.0, 1.0, 0.0, 0.0),
    vec4(0.0, 0.0, 1.0, 0.0),
    vec4(-eye.x, -eye.y, -eye.z, 1.0)
  );

  vec3 dir = normalize(at - eye);   
  mat4 R = mat4(
    vec4(X.x, Y.x, dir.x, 0.0),
    vec4(X.y, Y.y, dir.y, 0.0),
    vec4(X.z, Y.z, dir.z, 0.0),
    vec4(0,   0,   0,     1.0)
  );

  mat4 viewMatrix =  T * R;

  return viewMatrix;
}

mat4 g_vM;

void main() {
  // Get screen coordinates 
  // Origin at center of screen, range is [-1, 1]
  vec2 uv = 10.0*(gl_FragCoord.xy - 0.5*iRes.xy) / iRes.y;
  
  // get lookat matrix
  vec3 eye = vec3(5, 5, 5);
  vec3 at = vec3(0, 0, 0);
  vec3 up = vec3(0, 1, 0);
  g_vM = LookAt(eye, at, up);

  // set up ray 
  vec3 ro = vec3(6, 6, 6);     
  
  // pt on the screen positioned at the XY plane at Z=0
  //vec3 ps = vec3(uv.x, 0.0, uv.y);  
  //ps = (g_vM * vec4(ps.xyz, 1.0)).xyz;
  
  vec3 ps = vec3(uv.xy, 0.0);  
  ps = (g_vM * vec4(ps.xyz, 1.0)).xyz;

  vec3 rd = normalize(ps - ro); // ray direction 


  float d = rayMarch(ro, rd);
  vec3 p = ro + d * rd;
  
  float dif = getLight(p);

  // enable to debug normals
  // vec3 col = getNormal(p);
  float amb = 0.8;
  vec3 col = dif * vec3(1, 1, 1);
  fragColor = vec4(col, 1.0);

}
