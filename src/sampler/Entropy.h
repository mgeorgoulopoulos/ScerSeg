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
Probably not the best name was chosen here. The following class will take
samples from linear genome (which is a list of integers - each representing one
of the possible states). It also defines functions to comput Shannon entropy for
these fragments and also for random sets of integers picked randomly from the
whole genome.
*/

#include <QVector>
#include <random>

class EntropySampler {
  public:
	EntropySampler(const QVector<int> &pool, int windowSize = 25)
		: pool(pool), cloudDistribution(0, pool.size() - 1),
		  poolDistribution(0, pool.size() - windowSize - 1) {
		window.resize(windowSize);
		if (!pool.empty()) {
			int maxValue = *std::max_element(pool.begin(), pool.end());
			counts.resize(maxValue + 1);
			counts.fill(0);
		}
	}

	// Entropy of a random continuous slice of 'windowSize'
	double sampleSlice() const {
		const int start = poolDistribution(generator);
		return sampleSlice(start);
	}

	// Entropy of a given continuous slice of 'windowSize'
	double sampleSlice(int start) const {
		if (start >= pool.size() - window.size())
			throw(QString("Cannot sample from %1: not enough pool (size: %2, "
						  "window size: %3)")
					  .arg(start)
					  .arg(pool.size())
					  .arg(window.size()));
		for (int i = 0; i < window.size(); i++) {
			window[i] = pool[start + i];
		}
		return entropy(window);
	}

	// Entropy of a random set of genes from the pool
	double sampleCloud() const {
		for (int i = 0; i < window.size(); i++) {
			window[i] = pool[cloudDistribution(generator)];
		}
		return entropy(window);
	}

	double entropy(const QVector<int> &sample) const {
		counts.fill(0);
		for (const int v : sample) {
			Q_ASSERT(v >= 0 && v < counts.size());
			counts[v]++;
		}

		// calculate Shannon entropy
		double e = 0.0;
		for (const int c : counts) {
			if (c <= 0) {
				// zero probability
				continue;
			}
			// probability of the symbol
			const double p = (double)c / (double)sample.size();

			// Shannon formula
			e -= p * log2(p);
		}

		return e;
	}

  private:
	mutable QVector<int> counts;
	mutable QVector<int> window;
	const QVector<int> &pool;
	mutable std::default_random_engine generator;
	mutable std::uniform_int_distribution<int> cloudDistribution;
	mutable std::uniform_int_distribution<int> poolDistribution;
};
