#version 440

uniform mat4 MVMatrix, PMatrix;
uniform mat3 NormalMatrix;
uniform bool DisplayTruth;
uniform vec3 min_vertex, delta_vertex_inverse;

in vec2 TexCoord2D0;
in vec3 MCvertex, MCnormal, MCvertex_truth, MCnormal_truth;

out vec2 TexCoord2D;
out vec3 TexCoord3D, ecPosition, tnorm, ecPosition_truth, tnorm_truth;

out vec3 norm, norm_truth;		// 测试用

void main()
{
	TexCoord2D = TexCoord2D0;
	//TexCoord3D = (MCvertex - min_vertex) * delta_vertex_inverse;

	if (DisplayTruth)
	{
		ecPosition = vec3(MVMatrix * vec4(MCvertex_truth, 1.0f));
		gl_Position = PMatrix * MVMatrix * vec4(MCvertex_truth, 1.0f);
		tnorm = normalize(NormalMatrix * MCnormal_truth);
	}
	else
	{
		ecPosition = vec3(MVMatrix * vec4(MCvertex, 1.0f));
		gl_Position = PMatrix * MVMatrix * vec4(MCvertex, 1.0f);
		tnorm = normalize(NormalMatrix * MCnormal);
	}
	vec3 delta;
	float theta = 3.1415926535f / 3;
	delta.x = MCvertex.x - min_vertex.x;
	delta.y = MCvertex.y - min_vertex.y;
	delta.z = MCvertex.z - min_vertex.z;
	TexCoord3D.x = (cos(theta) * delta.x + sin(theta) * delta.y) / 0.3;
	TexCoord3D.y = (-sin(theta) * delta.x + cos(theta) * delta.y) / 0.3;
	TexCoord3D.z = delta.z / 0.3;
	ecPosition_truth = vec3(MVMatrix * vec4(MCvertex_truth, 1.0f));
	tnorm_truth = normalize(NormalMatrix * MCnormal_truth);

	norm = MCnormal;
	norm_truth = MCnormal_truth;
}
