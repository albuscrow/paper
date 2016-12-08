#version 440

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec2 te_g_texCoord2D[3];
in vec3 te_g_texCoord3D[3], te_g_ecPosition[3], te_g_tnorm[3], te_g_patchDistance[3];

//in float tePrimitive[3];
//in vec3 tePosition[3];
//in vec3 tePatchDistance[3];
//out vec3 gFacetNormal;
//out vec3 gPatchDistance;
//out vec3 gTriDistance;
//out float gPrimitive;

out vec2 TexCoord2D;
out vec3 TexCoord3D, ecPosition, tnorm, triDistance, patchDistance; 

void main()
{
	TexCoord2D = te_g_texCoord2D[0];
	TexCoord3D = te_g_texCoord3D[0];
	ecPosition = te_g_ecPosition[0];
	tnorm = te_g_tnorm[0];
    patchDistance = te_g_patchDistance[0];
    triDistance = vec3(1, 0, 0);
    gl_Position = gl_in[0].gl_Position; EmitVertex();

	TexCoord2D = te_g_texCoord2D[1];
	TexCoord3D = te_g_texCoord3D[1];
	ecPosition = te_g_ecPosition[1];
	tnorm = te_g_tnorm[1];
    patchDistance = te_g_patchDistance[1];
    triDistance = vec3(0, 1, 0);
    gl_Position = gl_in[1].gl_Position; EmitVertex();

	TexCoord2D = te_g_texCoord2D[2];
	TexCoord3D = te_g_texCoord3D[2];
	ecPosition = te_g_ecPosition[2];
	tnorm = te_g_tnorm[2];
    patchDistance = te_g_patchDistance[2];
    triDistance = vec3(0, 0, 1);
    gl_Position = gl_in[2].gl_Position; EmitVertex();

    EndPrimitive();
}
