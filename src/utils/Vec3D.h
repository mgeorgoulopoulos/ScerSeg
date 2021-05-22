/*
Copyright 2021 Michael Georgoulopoulos

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
Vector class that represents points in 3D space. Initially I used GLM, but I
decided it was not worth the dependency.
*/

#ifndef _VEC_3D_H_
#define _VEC_3D_H_

#include <cmath>

class Vec3D {
  public:
	Vec3D() {}
	Vec3D(double x, double y, double z) : x(x), y(y), z(z) {}

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	// Commonly used math operators
	Vec3D operator+(const Vec3D &other) const {
		return Vec3D(x + other.x, y + other.y, z + other.z);
	}

	void operator+=(const Vec3D &other) { *this = *this + other; }

	Vec3D operator-(const Vec3D &other) const {
		return Vec3D(x - other.x, y - other.y, z - other.z);
	}

	void operator-=(const Vec3D &other) { *this = *this - other; }

	Vec3D operator*(double multiplier) const {
		return Vec3D(x * multiplier, y * multiplier, z * multiplier);
	}

	void operator*=(double multiplier) { *this = *this * multiplier; }

	Vec3D operator/(double divisor) const {
		const double multiplier = 1.0 / divisor;
		return *this * multiplier;
	}

	void operator/=(double divisor) { *this = *this / divisor; }

	// Length (magnitude) of the vector
	double length() const { return sqrt(x * x + y * y + z * z); }

	// Euclidean distance
	static double distance(const Vec3D &a, const Vec3D &b) {
		const Vec3D displacement = b - a;
		return displacement.length();
	}

	// Linear interpolation of scalars
	static double mix(double a, double b, double t) {
		const double omt = 1.0 - t;
		return a * omt + b * t;
	}

	static Vec3D mix(const Vec3D &a, const Vec3D &b, double t) {
		return Vec3D(mix(a.x, b.x, t), mix(a.y, b.y, t), mix(a.z, b.z, t));
	}
};

#endif // _VEC_3D_H_
