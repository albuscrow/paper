#version 440

uniform sampler2D TexSamplerColormap, TexSampler2D;
uniform sampler3D TexSampler3D;
uniform samplerCube TexSamplerCube;
uniform bool UseEnvMap, ErrorTexture, LocalViewer;
uniform int TexCase, ErrorDisplayMode;
uniform mat3 CubeMapRotInvMat;

in vec2 TexCoord2D;
in vec3 TexCoord3D, ecPosition, tnorm, ecPosition_truth, tnorm_truth;

in vec3 norm, norm_truth;		// 测试用

out vec4 FragColor;

struct PointLightSource
{
	vec4 position, ambient, diffuse, specular;
	float constantAttenuation, linearAttenuation, quadraticAttenuation;
};
uniform PointLightSource PointLightSource0;

struct Material
{
	vec4 ka, kd, ks;
};
uniform Material Mtl;

void PointLight(const vec3 eye, const vec3 ecPosition3, const vec3 normal,
				inout vec4 ambient, inout vec4 diffuse, inout vec4 specular)
{
	vec3 VP = vec3(PointLightSource0.position) - ecPosition3;	// 片元到光源的向量
	float d = length(VP);										// 片元和光源的距离
	VP = normalize(VP);
	float attenuation = 1.0 / (PointLightSource0.constantAttenuation +
						 PointLightSource0.linearAttenuation * d +
						 PointLightSource0.quadraticAttenuation * d * d);// 衰减因子
// 方法1:
	float nDotHV = max(0.0, dot(reflect(-VP, normal), eye));
// 方法2:
	//halfVector = normalize(VP + eye);		// direction of maximum highlights
	//nDotHV = max(0.0, dot(normal, halfVector));

	float nDotVP = max(0.0, dot(normal, VP));	// normal.light
	float pf = 0.0;								// power factor
	if (nDotVP > 0.0)
		pf = pow(nDotHV, 10.0);

	ambient += PointLightSource0.ambient * attenuation;
	diffuse += PointLightSource0.diffuse * nDotVP * attenuation;
	specular += PointLightSource0.specular * pf * attenuation;
}

vec4 color_map_vertex(vec3 num, vec3 truth, float range)
{
	float len = distance(num, truth) / range;

	if (ErrorTexture)
		return texture(TexSamplerColormap, vec2(len, 0));
	else
	{
		vec4 color_0  = vec4(0.0, 0.0, 1.0, 1.0);
		vec4 color_05 = vec4(0.0, 1.0, 0.0, 1.0);
		vec4 color_1  = vec4(1.0, 0.0, 0.0, 1.0);
		if (len < 0.5)
			return mix(color_0, color_05, len * 2);
		else
			return mix(color_05, color_1, (len - 0.5) * 2);
	}
}

vec4 color_map_normal(vec3 norm0, vec3 norm1, float range)
{
	/*
	 * 下面两种误差计算方法的结果理论上相同，都是两个法向的夹角。
	 * 第一种不精确，即使norm0和norm1完全相同，得到的结果仍然不为0，
	 * 可能是因为点积时求各分量的乘积造成的；第二种较好
	 * 具体分析可以参看皮质笔记本
	 */
	//float arc_diff = acos(clamp(dot(norm0, norm1), -1.0, 1.0)) / range;
	//float arc_diff = 2 * asin(clamp(distance(norm0, norm1) * 0.5, 0.0, 1.0)) / range;

	// 理论上不需要clamp，因为这两个单位法向的距离的一半的范围本身就是[0.0, 1.0]
	// 2 * asin(distance(norm0, norm1) * 0.5)的范围是[0.0, PI]
	float arc_diff = 2 * asin(distance(norm0, norm1) * 0.5) / range;

	if (ErrorTexture)
		return texture(TexSamplerColormap, vec2(arc_diff, 0));
	else
	{
		vec4 color_0  = vec4(0.0, 0.0, 1.0, 1.0);
		vec4 color_05 = vec4(0.0, 1.0, 0.0, 1.0);
		vec4 color_1  = vec4(1.0, 0.0, 0.0, 1.0);
		if (arc_diff < 0.5)
			return mix(color_0, color_05, arc_diff * 2);
		else
			return mix(color_05, color_1, (arc_diff - 0.5) * 2);
	}
}

void main()
{
	if (ErrorDisplayMode == 1)			// 1表示显示顶点误差
		FragColor = color_map_vertex(ecPosition, ecPosition_truth, 0.02);
		//FragColor = color_map_vertex(ecPosition, ecPosition_truth, 0.0007);
	else if (ErrorDisplayMode == 2)		// 2表示显示法向误差
		//FragColor = vec4(tnorm - tnorm_truth, 1.0);
		//FragColor = vec4(norm - norm_truth, 1.0);
		FragColor = color_map_normal(tnorm, tnorm_truth, 3.14159265358979 / 70);
		//FragColor = color_map_normal(tnorm, tnorm_truth, 3.14159265358979 / 5);
		//FragColor = color_map_normal(tnorm, tnorm_truth, 3.14159265358979 / 4);
	else
	{
		vec3 eye = vec3(0.0, 0.0, 1.0);
		if (LocalViewer)
			eye = -normalize(ecPosition);
		vec4 amb = vec4(0.0), diff = vec4(0.0), spec = vec4(0.0);
		PointLight(eye, ecPosition, normalize(tnorm), amb, diff, spec);

		FragColor = amb * Mtl.ka + diff * Mtl.kd;
		vec4 SecondaryColor = spec * Mtl.ks;
		FragColor = clamp(FragColor, 0.0, 1.0);

		if (TexCase == 2)
		{
			/*FragColor = texture(TexSampler2D, TexCoord2D);*/
			FragColor *= texture(TexSampler2D, TexCoord2D);
		}
		else if (TexCase == 3)
		{
			FragColor *= texture(TexSampler3D, TexCoord3D);
		}
		if (UseEnvMap)
		{
			vec3 ReflectDir = CubeMapRotInvMat * reflect(ecPosition, tnorm);
			ReflectDir.z = -ReflectDir.z;
			FragColor = texture(TexSamplerCube, ReflectDir);
		}
		FragColor += SecondaryColor;
		FragColor = clamp(FragColor, 0.0, 1.0);

		/*FragColor = vec4(0.0, ecPosition.x, 1.0 - ecPosition.x, 1.0);*/
		/*FragColor = vec4(ecPosition.x, 1.0 - ecPosition.x, 0.0, 1.0);*/
		/*FragColor = SecondaryColor;*/
		/*FragColor = vec4(norm, 1.0);*/
		/*FragColor = vec4(pow(norm.x, 10.0), pow(norm.y, 10.0), pow(norm.z, 10.0), 1.0);*/
		/*FragColor = vec4(pow(tnorm.x, 10.0), pow(tnorm.y, 10.0), pow(tnorm.z, 10.0), 1.0);*/
		/*FragColor = vec4(abs(eye), 1.0);*/
		/*FragColor = abs(PointLightSource0.position) + 0.5;*/
		/*FragColor = vec4(normalize(norm), 1.0);*/
	}
}
