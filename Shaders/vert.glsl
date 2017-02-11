#version 450 core

layout (location=0) in vec3 position;
layout (location=1) in vec3 normal;
layout (location=2) uniform float t;
layout (location=3) uniform uint coloring;
layout (location=4) uniform float lowest;
layout (location=5) uniform float highest;
layout (location=6) uniform float dist;

uniform vec3 lightIntensities = vec3(0.8)/*(0.5, 0.4195, 0.3)*/;
uniform vec3 ambient = vec3(0.8);

out VS_OUT
{   vec4 color;
    vec3 normal;
    vec3 pos;
    vec3 lightPos;
} vs_out;

uniform mat4 pMat;
uniform mat4 vMat;
uniform mat4 mMat;

void main(void)
{
    /* backwards interpolation */
    float k = (position.y-lowest)/(highest-lowest);
    gl_Position = pMat * vMat * mMat * vec4(position,1.0);
    float ns = clamp(noise1(position.y), 0.2, 0.6);
    if (normal[0] == 2)
        vs_out.color = vec4(0.34901960784313724, 0.8431372549019608, 0.0, 1.0);
    else if (coloring == 0)
        vs_out.color = k <= 0.3 ? vec4(0.2, 0.2, 1.0, 1.0) : vec4(vec3(0.5*k), 1.0); 
    else if (coloring == 1)
        vs_out.color = vec4(1.0-clamp(.6/(abs(position.y)+0.01), 0.0, 1.0), k, clamp(.4/(abs(position.y)+0.01), 0.0, 1.0), 1.0);
    else if (coloring == 2)
        vs_out.color = vec4(abs(1./exp(position.y)), abs(cos(5*position.y)), abs(sin(5*position.y)), 1.0);
    else if (coloring == 3)
        vs_out.color = vec4(ns)+mix(vec4(0.0, 0.0, 0.8, 1.0), vec4(0.5, 1.0, 0.3, 1.0), k);
    else if (coloring == 4)
        vs_out.color = vec4(ns)+mix(vec4(1.0), vec4(0.1, 0.3, 1.0, 0.5), k);
    vs_out.normal = normalize(mat3(vMat*mMat)*normal);
    vs_out.pos = mat3(vMat*mMat) * position; 
    vs_out.lightPos = gl_Position.xyz;

    vec3 l = normalize(vs_out.lightPos-vs_out.pos);
    float brightness = clamp(dot(vs_out.normal, l), 0.0, 1.0);
    vs_out.color = vec4((ambient + brightness * lightIntensities) * vs_out.color.rgb, vs_out.color.a);
}
