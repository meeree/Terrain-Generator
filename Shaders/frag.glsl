#version 450 core

out vec4 color;

in VS_OUT
{   vec4 color;
    vec3 normal;
    vec3 pos;
    vec3 lightPos;
} fs_in;

void main(void)
{
    color = fs_in.color;
}
