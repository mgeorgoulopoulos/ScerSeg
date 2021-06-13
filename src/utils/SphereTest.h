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
This module contains templated functions that facilitate the sphere test,
described in detail in the accompanying documentation. All templates get a Gene
argument which corresponds to the gene class of the particular instance of the
test. This is an effort to reuse parts of the algorithm in multiple programs.
The test itself is not templated and is a (copied) large function in each
specific program. We at least deflate the large function using these templates
here. We dont include the large function here so the various programs can have
some room to recombine parts of the procedure.
*/

#ifndef _SPHERE_TEST_H_
#define _SPHERE_TEST_H_

// One work unit is one successfully-sampled sphere together with its results
// (genes and p-values). We use a separate struct for this to facilitate
// multithreading. Each sphere does its p-value calculation in its own thread.
// Each thread has its own private data to work with.
template <class Gene> struct WorkUnit {
	QVector<const Gene *> genesInSphere;
	double pValue = 1.0;
	Sampler::RandomGeneSampler<Gene> randomSampler;

	// How many times we got a more extreme result by sheer luck?
	int chanceWinCount = 0;

	// Buffer for sampling random sets
	QVector<const Gene *> randomGenes;

	// Keep these 2 for further processing - random is the last of the random tests
	double statisticInSphere = 0.0;
	double statisticInRandom = 0.0;

	// P-value rank and Benjamini-Hochberg-adjusted version of it.
	int rank = 0;
	double adjustedPValue = 1.0;

	WorkUnit() {}
	WorkUnit(const WorkUnit &other) : randomSampler(other.randomSampler) {
		genesInSphere = other.genesInSphere;
		pValue = other.pValue;
		chanceWinCount = other.chanceWinCount;
		randomGenes = other.randomGenes;
		rank = other.rank;
		adjustedPValue = other.adjustedPValue;
	}

	WorkUnit(const QVector<Gene> &pool) : randomSampler(pool) {}

	void calculatePValue(int randomSampleCount) {
		// Calculate metric in the sphere-sample
		statisticInSphere = sphereTestStatistic(genesInSphere);

		// Random samples
		for (int r = 0; r < randomSampleCount; r++) {
			randomSampler.sample(genesInSphere.size(), &randomGenes);
			statisticInRandom = sphereTestStatistic(randomGenes);

			// Is random more extreme than sphere?
			if (Gene::randomIsMoreExtreme(statisticInRandom, statisticInSphere))
				chanceWinCount++;
		} // end for (N random samples)

		pValue = Gene::calculatePValue(chanceWinCount, randomSampleCount);

		// Give benefit of the doubt to chance: replace zero pValues with the
		// smallest we can safely say
		pValue = std::max(pValue, 1.0 / (double)randomSampleCount);
	}
};

// Creates a number of randomized WorkUnits.
template <class Gene>
QVector<WorkUnit<Gene>> createWorkUnits(double sphereRadius,
										int minimumGeneCount,
										const QVector<Gene> &genes, int count,
										int *averageGenesInASphere) {
	QVector<WorkUnit<Gene>> result;
	result.reserve(count);

	const Sampler::SphereGeneSampler<Gene> sphereSampler(genes, boxMinimum,
														 boxMaximum);

	*averageGenesInASphere = 0;

	while (result.size() < count) {
		WorkUnit<Gene> workUnit(genes);

		workUnit.genesInSphere = sphereSampler.sample(sphereRadius);
		if (workUnit.genesInSphere.size() < minimumGeneCount) {
			// Reject sample - we don't want spheres in mostly empty space
			continue;
		}

		// Successful sample
		result.push_back(workUnit);
		*averageGenesInASphere += workUnit.genesInSphere.size();
	}

	*averageGenesInASphere /= count;

	return result;
}

// Adjusts p-values and reorders all work units.
template <class Gene> void benjamini(QVector<WorkUnit<Gene>> &workUnits) {
	// Sort p-values from smaller to larger
	auto pValueCompare = [](const WorkUnit<Gene> &a, const WorkUnit<Gene> &b) {
		return a.pValue < b.pValue;
	};
	std::sort(workUnits.begin(), workUnits.end(), pValueCompare);

	// Apply Benjamini-Hochberg correction
	double previousPValue = workUnits.back().pValue;
	for (int i = workUnits.size() - 1; i >= 0; i--) {
		WorkUnit<Gene> &workUnit = workUnits[i];
		workUnit.rank = i + 1;
		workUnit.adjustedPValue =
			std::min(previousPValue, workUnit.pValue * workUnits.size() /
										 (double)workUnit.rank);
		previousPValue = workUnit.adjustedPValue;
	}
}

// Given a list of work units, it combines the overlapping spheres into clusters
template <class Gene>
QVector<QSet<QString>>
clusterByGeneOverlap(const QVector<WorkUnit<Gene>> &workUnits,
					 double overlapThreshold, double *maximumOverlapRatio) {
	using GeneList = QSet<QString>;
	QVector<GeneList> result;
	for (const WorkUnit<Gene> &workUnit : workUnits) {
		GeneList genes;
		for (const Gene *gene : workUnit.genesInSphere) {
			genes.insert(gene->name);
		}
		result.push_back(genes);
	}

	// Now do hierarchical clustering.
	while (true) {
		// Find max-overlap pair of clusters
		*maximumOverlapRatio = 0.0;
		using GenePair = QPair<int, int>;
		GenePair bestGenePair(0, 1);
		GeneList mergedCluster;
		for (int i = 0; i < result.size(); i++) {
			for (int j = i + 1; j < result.size(); j++) {
				const int minSize =
					std::min(result[i].size(), result[j].size());
				GeneList intersection = result[i];
				intersection.intersect(result[j]);
				const int overlapCount = intersection.size();
				const double overlapRatio =
					(double)overlapCount / (double)minSize;
				if (overlapRatio > *maximumOverlapRatio) {
					*maximumOverlapRatio = overlapRatio;
					mergedCluster = result[i];
					mergedCluster.unite(result[j]);
					bestGenePair = GenePair(i, j);
				}
			}
		} // end for (all pairs)

		// When we reach the point where the best-overlapping clusters may be
		// considered distinct, we are done.
		if (*maximumOverlapRatio < overlapThreshold) {
			break;
		}

		// Not done yet: apply the merge and continue.
		result[bestGenePair.second] = result.back();
		result.pop_back();
		result[bestGenePair.first] =
			result.back(); // we can do this always: first is less than second
						   // by design.
		result.pop_back();
		result.push_back(mergedCluster);
	} // end hierarchical clustering

	return result;
}

#endif // _SPHERE_TEST_H_
