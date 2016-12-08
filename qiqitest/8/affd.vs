#version 440

uniform mat4 MVMatrix, PMatrix;
uniform mat3 NormalMatrix;

in vec2 TexCoord2D0;
in vec3 TexCoord3D0, MCvertex, MCnormal; 

out vec2 TexCoord2D;
out vec3 TexCoord3D, ecPosition, tnorm; 

void main()
{
	TexCoord2D = TexCoord2D0;
	TexCoord3D = TexCoord3D0;

	ecPosition = vec3(MVMatrix * vec4(MCvertex, 1.0f));
	gl_Position = PMatrix * MVMatrix * vec4(MCvertex, 1.0f);
	tnorm = normalize(NormalMatrix * MCnormal);

	/************** 以下为旋转三位纹理的代码，如无必要可以注释掉 ***************/
	//float theta = 3.1415926535f * 0.333333333;
	float theta = 3.1415926535f * 0.5;
	//float theta = 3.1415926535f / 3;
	float factor = 1;

	// 绕x轴旋转
	TexCoord3D.x = TexCoord3D0.x * factor;
	TexCoord3D.y = (cos(theta) * TexCoord3D0.y - sin(theta) * TexCoord3D0.z) * factor;
	TexCoord3D.z = (sin(theta) * TexCoord3D0.y + cos(theta) * TexCoord3D0.z) * factor;

	// 绕y轴旋转
	//TexCoord3D.x = (cos(theta) * TexCoord3D0.x + sin(theta) * TexCoord3D0.z) * factor;
	//TexCoord3D.y = TexCoord3D0.y * factor;
	//TexCoord3D.z = (-sin(theta) * TexCoord3D0.x + cos(theta) * TexCoord3D0.z) * factor;

	// 绕z轴旋转
	//TexCoord3D.x = (cos(theta) * TexCoord3D0.x - sin(theta) * TexCoord3D0.y) * factor;
	//TexCoord3D.y = (sin(theta) * TexCoord3D0.x + cos(theta) * TexCoord3D0.y) * factor;
	//TexCoord3D.z = TexCoord3D0.z * factor;
}
