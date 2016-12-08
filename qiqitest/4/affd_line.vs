#version 440

in vec2 TexCoord2D0;
in vec3 v_ctrlPoint, v_ctrlPointN; 

out vec2 v_tc_texCoord2D;
out vec3 v_tc_ctrlPoint, v_tc_ctrlPointN;

void main()
{
	//TexCoord2D = TexCoord2D0;
	//TexCoord3D = TexCoord3D0;

	//ecPosition = vec3(MVMatrix * vec4(MCvertex, 1.0f));
	//gl_Position = PMatrix * MVMatrix * vec4(MCvertex, 1.0f);
	//tnorm = normalize(NormalMatrix * MCnormal);

	v_tc_texCoord2D = TexCoord2D0;
	v_tc_ctrlPoint = v_ctrlPoint;
	v_tc_ctrlPointN = v_ctrlPointN;
}
