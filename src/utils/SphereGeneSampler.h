/*
Copyright 2021 Michael Georgoulopoulos

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files(the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and / or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
This module defines templated class SphereSampler. SphereSampler was a piece of
code I copied and pasted a lot so it probably deserves its own module. The
reason for not doing so from the start was that "Gene" meant very different
things to different programs. Some program might have histone modification
profiles included in the Gene struct, some other conservation, transcription
factor motifs and so on. The natural way of handling this is a C++ template
which requires the common denominator of all structs: having a name and a 3D
position so we can sample them.
*/

#ifndef _SPHERE_SAMPLER_H_
#define _SPHERE_SAMPLER_H_

#include <QVector>
#include <random>

#include "Vec3D.h"

namespace Sampler {

template <class Gene> class SphereGeneSampler {
  public:
	SphereGeneSampler(const QVector<Gene> &pool, double boxMinimum,
					  double boxMaximum)
		: pool(pool), distribution(boxMinimum, boxMaximum),
		  boxMinimum(boxMinimum), boxMaximum(boxMaximum) {}

	// Samples genes from the pool and returns pointers to them, using the
	// specified radius and a random center. The center x,y,z coordinates are
	// taken using a uniform distribution in the user-provided box.
	QVector<const Gene *> sample(double radius) const {
		QVector<const Gene *> result;

		Vec3D center;
		center.x = distribution(generator);
		center.y = distribution(generator);
		center.z = distribution(generator);

		for (const Gene &gene : pool) {
			if (Vec3D::distance(center, gene.position) <= radius) {
				result.push_back(&gene);
			}
		}

		return result;
	}

  private:
	const QVector<Gene> &pool;
	const double boxMinimum;
	const double boxMaximum;
	mutable std::default_random_engine generator;
	mutable std::uniform_real_distribution<double> distribution;
};

} // end namespace Sampler
#endif // _SPHERE_SAMPLER_H_
