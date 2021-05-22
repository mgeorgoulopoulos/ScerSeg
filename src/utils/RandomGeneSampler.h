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
Random gene sampler. Template similar to SphereGeneSampler. Provided for
symmetry. Returns fully random Gene sets with replacement from the pool
provided.
*/

#ifndef _RANDOM_GENE_SAMPLER_H_
#define _RANDOM_GENE_SAMPLER_H_

#include <QVector>
#include <random>

namespace Sampler {

template <class Gene> class RandomGeneSampler {
  public:
	RandomGeneSampler() {}
	RandomGeneSampler(const QVector<Gene> &pool)
		: pool(&pool), generator(std::random_device()()),
		  distribution(0, pool.size() - 1) {}
	RandomGeneSampler(const RandomGeneSampler &other)
		: pool(other.pool), generator(other.generator),
		  distribution(other.distribution) {}

	void sample(int count, QVector<const Gene *> *genes) const {
		genes->resize(0);

		if (pool->isEmpty())
			return;

		for (int i = 0; i < count; i++) {
			genes->push_back(&(*pool)[distribution(generator)]);
		}
	}

	int random() const { return distribution(generator); }

  private:
	const QVector<Gene> *pool = nullptr;
	mutable std::default_random_engine generator;
	mutable std::uniform_int_distribution<int> distribution;
};

} // end namespace Sampler

#endif // _RANDOM_GENE_SAMPLER_H_
