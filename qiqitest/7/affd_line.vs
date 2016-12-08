#version 440

uniform mat4 MVMatrix, PMatrix;
uniform mat3 NormalMatrix;
uniform vec3 min_vertex, delta_vertex_inverse;

in vec2 TexCoord2D0;
in vec3 MCvertex, MCnormal, Bary, OriBary;

out vec2 TexCoord2D;
out vec3 TexCoord3D, ecPosition, tnorm, bary, oriBary;

void main()
{
	TexCoord2D = TexCoord2D0;

	//TexCoord3D = (MCvertex - min_vertex) * delta_vertex_inverse;

	//TexCoord3D = (MCvertex - min_vertex) / 0.3;

	vec3 delta;
	float theta = 3.1415926535f / 3;
	delta.x = MCvertex.x - min_vertex.x;
	delta.y = MCvertex.y - min_vertex.y;
	delta.z = MCvertex.z - min_vertex.z;
	TexCoord3D.x = (cos(theta) * delta.x + sin(theta) * delta.y) / 0.3;
	TexCoord3D.y = (-sin(theta) * delta.x + cos(theta) * delta.y) / 0.3;
	TexCoord3D.z = delta.z / 0.3;

	//TexCoord3D.x = (MCvertex.y - min_vertex.y) / 0.3;
	//TexCoord3D.y = (MCvertex.z - min_vertex.z) / 0.3;
	//TexCoord3D.z = (MCvertex.x - min_vertex.x) / 0.3;

	bary = Bary;
	oriBary = OriBary;

	ecPosition = vec3(MVMatrix * vec4(MCvertex, 1.0f));
	gl_Position = PMatrix * MVMatrix * vec4(MCvertex, 1.0f);
	tnorm = normalize(NormalMatrix * MCnormal);
}
