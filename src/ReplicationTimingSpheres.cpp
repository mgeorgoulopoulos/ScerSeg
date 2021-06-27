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
This program implements the sphere sampling algorithm discussed in accompanying
documentation. We use replication timing as input.
*/

// Settings
namespace {

// Description of the test. This will appear in STDOUT
const char *statisticDescription = "Standard deviation of replication timing in sphere.";

// Table name in the DB, where the resulting clusters will be persisted.
const char *tableName = "ReplicationTImingFields";

// Radius of the sampling sphere
const double sphereRadius = 15.0;

// Threshold to disregard mostly empty spheres
const int minimumGeneCount = 50;

// Minimum/maximum x,y,z dimension of the gene box. This box is where random
// spheres are picked from.
const double boxMinimum = 0.0;
const double boxMaximum = 210.0;

// Pay attention to this: it is how many sphere samples we get and
// subsequently: how many random samples of the same size for each
// sphere. So processing time is O(n^2) to this.
const int sampleCount = 10000;

// Filter sphere samples by adjusted p-value. I propose to run this program
// twice: on first run p-values can be examined (they are written to text file).
// Subsequently, you can set this to a sane value, so that only significant
// spheres are taken into consideration. Typical values range from 1% to 5%.
const double pAdjThreshold = 0.05;

// Used for hierarchical clustering of overlapping spheres into the
// final result. Less than this overlap defines two distinct clusters.
// To avoid fairly separate clusters become a huge blob for just
// sharing one or a few genes.
const double overlapThreshold = 0.05;

} // namespace

#include "utils/RandomGeneSampler.h"
#include "utils/SaveClustersToDB.h"
#include "utils/SphereGeneSampler.h"
#include "utils/SphereTest.h"
#include "utils/Vec3D.h"

#include <QElapsedTimer>
#include <QFile>
#include <QSet>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QSqlRecord>
#include <QString>
#include <QTextStream>
#include <QVector>

#include <algorithm>
#include <bitset>
#include <random>

#define TF_COUNT 102

namespace {

// Define your Gene structure here. It is required that Gene at least contains:
//		* QString name
//		* Vec3D position
//		* p-value helper functions. See below. These basically define what
//			consists an "extreme" outcome.
struct Gene {
	QString name;
	Vec3D position;
	double replicationTiming = 0.0;
	int orderInGenome;

	// Implement this function to define what is considered more extreme
	static bool randomIsMoreExtreme(double randomStatistic,
									double testStatistic) {
		// More varied than random seems to be the rare thing here.
		return randomStatistic >= testStatistic;
	}

	// Implement this to calculate p-value. Single- and Double-tailed tests may
	// be implemented this way.
	static double calculatePValue(int randomMoreExtremeOccurrences,
								  int totalRandomSamples) {
		// Calculate 1-tailed p-value
		double pValue =
			(double)randomMoreExtremeOccurrences / (double)totalRandomSamples;

		// Use this to convert to 2-tailed test
		// pValue = 2.0 * std::min(pValue, 1.0 - pValue);

		return pValue;
	}

	// Define this to control when a sphere-sample is acceptable
	static bool acceptSample(const QVector<const Gene *> &genes) {
		if (genes.size() < minimumGeneCount)
			return false;

		// Replication timing is a smooth signal. Therefore we want spheres that
		// conver more than one chromosome, or more than one locations of the
		// same chromosome. We encode this as a gene index jump of at least 100
		// genes.
		QVector<const Gene *> sorted = genes;
		std::sort(sorted.begin(), sorted.end(),
				  [](const Gene *a, const Gene *b) {
					  return a->orderInGenome < b->orderInGenome;
				  });
		const int minimumIndexSpaceJump = 100;
		for (int i = 1; i < sorted.size(); i++) {
			if (sorted[i]->orderInGenome - sorted[i - 1]->orderInGenome >=
				minimumIndexSpaceJump)
				return true;
		}

		return false;
	}
};

// Load set of genes from DB
QVector<Gene> loadGenes(QSqlDatabase &db, const QString &whereClause = "") {
	QVector<Gene> result;

	const QString sql = QString(
		"SELECT l.Gene, x, y, z, ReplicationTiming FROM Loci l JOIN "
		"ReplicationTiming r ON l.Gene = r.Gene ORDER BY Chromosome, Start");
	QSqlQuery query(sql, db);
	while (query.next()) {
		Gene gene;
		gene.name = query.value(0).toString();
		gene.position.x = query.value(1).toDouble();
		gene.position.y = query.value(2).toDouble();
		gene.position.z = query.value(3).toDouble();
		gene.replicationTiming = query.value(4).toDouble();
		gene.orderInGenome = result.size();

		result.push_back(gene);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

// This implements the sphere test statistic for our particular test.
double sphereTestStatistic(const QVector<const Gene *> &genes) {
	if (genes.isEmpty())
		throw(QString("Empty gene list provided!"));

	// Calculate average
	double average = 0.0;
	for (const Gene *gene : genes) {
		average += gene->replicationTiming;
	}
	average /= (double)genes.size();

	double variance = 0.0;
	for (const Gene *gene : genes) {
		double averageDiff = gene->replicationTiming - average;
		variance += averageDiff * averageDiff;
	}
	variance /= (double)genes.size();

	const double stdev = sqrt(variance);

	return stdev;
}

// This function performs the sphere test. An almost identical copy of this
// function exists in all similar sphere-test programs. This is a compromise
// between reusability and flexibility. Heavy parts of the procedure have been
// extracted to SphereTest.h.
void extractReplicationTimingFields(QSqlDatabase &db) {
	printf("Description of measured statistic:\n\t%s\n", statisticDescription);

	// Load genes
	const QVector<Gene> genes = loadGenes(db);
	printf("%d genes\n", genes.size());

	QElapsedTimer timer;
	timer.start();

	// Generate a number of sphere samples

	printf("Generating %d sphere samples... ", sampleCount);
	int averageGenesInASphere = 0;
	QVector<WorkUnit<Gene>> workUnits = createWorkUnits(
		sphereRadius, genes, sampleCount, &averageGenesInASphere);
	printf("Done.\n");
	printf("Average genes in a sphere: %d\n", averageGenesInASphere);

	// Then, for each of the samples, draw the same number of random samples of
	// the same gene count. Calculate the same metric and calculate a p-value.

	// Reuse random samples. We don't really need thousands of random samples
	// for each of the thousands of sphere samples. We need them for each *gene
	// count*.
	printf("Calculating statistic on %d random samples for all possible gene "
		   "set sizes...",
		   workUnits.size());
	QMap<int, QVector<double>> geneCountToRandomStatistics;
	for (WorkUnit<Gene> &workUnit : workUnits) {
		const int geneCount = workUnit.genesInSphere.size();
		if (geneCountToRandomStatistics.contains(geneCount))
			continue;
		printf("%d ", geneCount);

		geneCountToRandomStatistics[geneCount].reserve(sampleCount);
		for (int r = 0; r < sampleCount; r++) {
			workUnit.randomSampler.sample(geneCount, &workUnit.randomGenes);
			const double statisticInRandom =
				sphereTestStatistic(workUnit.randomGenes);
			geneCountToRandomStatistics[geneCount].push_back(statisticInRandom);
		}
	}
	printf("Done\n");

	printf("Calculating p-value for each of %d spheres... ", workUnits.size());
	for (int i = 0; i < workUnits.size(); i++) {
		WorkUnit<Gene> &workUnit = workUnits[i];

		// Calculate metric in the sphere-sample
		const double statisticInSphere =
			sphereTestStatistic(workUnit.genesInSphere);

		// Random samples
		const QVector<double> &statsOnRandomSample =
			geneCountToRandomStatistics[workUnit.genesInSphere.size()];
		for (const double statisticInRandom : statsOnRandomSample) {
			// Is random more extreme than sphere?
			if (Gene::randomIsMoreExtreme(statisticInRandom, statisticInSphere))
				workUnit.chanceWinCount++;
		} // end for (N random samples)

		workUnit.pValue =
			Gene::calculatePValue(workUnit.chanceWinCount, sampleCount);

		// Give benefit of the doubt to chance: replace zero pValues with the
		// smallest we can safely say
		workUnit.pValue = std::max(workUnit.pValue, 1.0 / (double)sampleCount);
	}
	printf("Done.\n");

	// Adjust p-values
	printf("Adjusting p-values using Benjamini-Hochberg method... ");
	benjamini(workUnits);
	printf("Done.\n");

	// Write p-values to file for later reference
	{
		const QString filename =
			QString("Results/pValues.%1.tsv").arg(tableName);
		printf("Writing p-values to file: %s\n", filename.toUtf8().data());
		QFile file(filename);
		if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
			throw QString("Failed to open file %1 for writing\n").arg(filename);
		QTextStream out(&file);
		for (const WorkUnit<Gene> &workUnit : workUnits) {
			out << workUnit.pValue << '\t' << workUnit.adjustedPValue << '\n';
		}
		file.close();
	}

	// Filter by adjusted p-value
	QVector<WorkUnit<Gene>> tmp;
	std::copy_if(workUnits.begin(), workUnits.end(), std::back_inserter(tmp),
				 [&](const WorkUnit<Gene> &x) {
					 return x.adjustedPValue <= pAdjThreshold;
				 });
	workUnits = tmp;
	printf("%d significant p-values (below %.05f)\n", workUnits.size(),
		   pAdjThreshold);

	// Get the significant genes
	QSet<QString> significantGenes;
	for (const WorkUnit<Gene> &workUnit : workUnits) {
		for (const Gene *gene : workUnit.genesInSphere) {
			significantGenes.insert(gene->name);
		}
	}
	printf("%d significant genes.\n", significantGenes.size());

	if (!workUnits.isEmpty()) {
		printf("Here is the best sphere sample: (p-value: %f)\n",
			   workUnits.front().pValue);
		for (const Gene *gene : workUnits.front().genesInSphere) {
			printf("%s ", gene->name.toUtf8().data());
		}
		printf("\n");
	} else {
		printf("No significant samples found");
	}

	// Convert work units to clusters of genes
	printf("\n");
	printf("Hierarchical clustering. Using threshold of %.02f%% overlap ratio "
		   "to consider clusters as distinct ... ",
		   overlapThreshold * 100.0);
	double maximumOverlapRatio = 0.0;
	QVector<QSet<QString>> clusters =
		clusterByGeneOverlap(workUnits, overlapThreshold, &maximumOverlapRatio);
	printf("Done.\n");
	printf("Stopping clustering with %d clusters, %.02f%% maximum gene "
		   "overlap.\n",
		   clusters.size(), maximumOverlapRatio * 100.0);

	// The returned clusters are fairly defined. Sometimes we get a different
	// number of clusters because this is a random sampling, but mostly we get
	// the same. Ordering the clusters by size.
	std::sort(clusters.begin(), clusters.end(),
			  [](const QSet<QString> &a, const QSet<QString> &b) {
				  return a.size() < b.size();
			  });

	// Cleanup duplicate genes - let the smaller cluster retain all its genes
	QSet<QString> genesUsed;
	for (QSet<QString> &cluster : clusters) {
		cluster.subtract(genesUsed);
		genesUsed.unite(cluster);
	}

	// Report clusters
	QMap<QString, Gene> nameToGene;
	for (const Gene &gene : genes) {
		nameToGene.insert(gene.name, gene);
	}
	for (int i = 0; i < clusters.size(); i++) {
		QVector<const Gene *> tmp;
		for (const QString &geneName : clusters[i]) {
			tmp.push_back(&nameToGene[geneName]);
		}
		const double metric = sphereTestStatistic(tmp);

		printf("\tCluster %d: %d genes\tMetric=%f\n", i + 1, clusters[i].size(),
			   metric);
	}

	// Calculate the metric for the entire population to give a hint on range
	{
		QVector<const Gene *> tmp;
		for (const Gene &gene : genes) {
			tmp.push_back(&gene);
		}
		const double metricOverGenome = sphereTestStatistic(tmp);
		printf("Metric calculated over the entire genome = %f\n",
			   metricOverGenome);
	}

	// Write clusters to database for further evaluation
	if (!clusters.isEmpty()) {
		printf("Writing clusters to database ... ");
		db::writeClusters(clusters, db, tableName);
		printf("Done.\n");
	} else {
		printf("No clusters found - not creating a table\n");
	}

	// Done. Report time it took
	printf("\nElapsed time: %.02f minutes.\n", timer.elapsed() / 60000.0);
}

} // end anonymous namespace

int main(int argc, char *argv[]) {
	const QString &filename = "Results/yeast.sqlite";

	if (!QFile::exists(filename)) {
		printf("No such file: %s\n", filename.toUtf8().data());
		return 0;
	}

	QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE");
	db.setDatabaseName(filename);
	if (!db.open()) {
		printf("Failed to open file: %s\n", filename.toUtf8().data());
	}

	try {
		extractReplicationTimingFields(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
