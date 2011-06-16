#pragma once

struct Color4 {

	Color4(GLubyte color = 255){
		m_v[0] = color;
		m_v[1] = color;
		m_v[2] = color;
		m_v[3] = 255;
	}

	Color4(GLfloat *v){
		m_v[0] = (GLubyte) (v[0] * 255);
		m_v[1] = (GLubyte) (v[1] * 255);
		m_v[2] = (GLubyte) (v[2] * 255);
		m_v[3] = 255;
	}

	Color4(GLubyte r, GLubyte g, GLubyte b, GLubyte a = 255){
		m_v[0] = r;
		m_v[1] = g;
		m_v[2] = b;
		m_v[3] = a;
	}

	static Color4 fromRGB(float r, float g, float b, float a = 1.0)
	{
		Color4 c;

		c.m_v[0] = r * 255;
		c.m_v[1] = g * 255;
		c.m_v[2] = b * 255;
		c.m_v[3] = a * 255;

		return c;
	}

	static Color4 random()
	{
		Color4 c;

		c.m_v[0] = rand() % 255;
		c.m_v[1] = rand() % 255;
		c.m_v[2] = rand() % 255;
		c.m_v[3] = 255;

		return c;
	}

	void set(float r, float g, float b){
		m_v[0] = (GLubyte) (r * 255);
		m_v[1] = (GLubyte) (g * 255);
		m_v[2] = (GLubyte) (b * 255);
		m_v[3] = 255;
	}

	inline void set(int r, int g, int b, int a = 255){
		m_v[0] = (GLubyte)r;
		m_v[1] = (GLubyte)g;
		m_v[2] = (GLubyte)b;
		m_v[3] = (GLubyte)a;
	}

	inline float r() const { return (float) m_v[0] / 255.0; }
	inline float g() const { return (float) m_v[1] / 255.0; }
	inline float b() const { return (float) m_v[2] / 255.0; }
	inline float a() const { return (float) m_v[3] / 255.0; }

	GLubyte operator[](int i){ 
		return m_v[i];
	}

	operator GLubyte*()
	{
		return m_v;
	}

	GLubyte m_v[4];
};
