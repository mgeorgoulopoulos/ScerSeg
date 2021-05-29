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
documentation. We use number of species having each gene as input.
*/

// Define / undefine this to switch between taxon and species count tests
#define TAXON_TEST

// Settings
namespace {

// Description of the test. This will appear in STDOUT
#ifdef TAXON_TEST
	const char *statisticDescription =
		"Maximum taxon absolute log enrichment.";
#else
	const char *statisticDescription =
	"Standard deviation in number of species having a gene.";
#endif

// Table name in the DB, where the resulting clusters will be persisted.
#ifdef TAXON_TEST
	const char *tableName = "TaxonFields";
#else
const char *tableName = "ConservationFields";
#endif

// Radius of the sampling sphere
const double sphereRadius = 10.0;

// Threshold to disregard mostly empty spheres
const int minimumGeneCount = 20;

// Minimum/maximum x,y,z dimension of the gene box. This box is where random
// spheres are picked from.
const double boxMinimum = 0.0;
const double boxMaximum = 210.0;

// Pay attention to this: it is how many sphere samples we get and
// subsequently: how many random samples of the same size for each
// sphere. So processing time is O(n^2) to this.
const int sampleCount = 20000;

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

#define HISTONE_COLUMN_COUNT 9

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

namespace {

// Define your Gene structure here. It is required that Gene at least contains:
//		* QString name
//		* Vec3D position
//		* p-value helper functions. See below. These basically define what
//			consists an "extreme" outcome.
struct Gene {
	QString name;
	Vec3D position;
	int speciesCount;
	QString taxon;

	// Implement this function to define what is considered more extreme
	static bool randomIsMoreExtreme(double randomStatistic,
									double testStatistic) {
		
#ifdef TAXON_TEST
		// Taxon enrichment: greater is considered extreme.
		return randomStatistic >= testStatistic;
#else
		// Standard deviation of species count: less is considered extreme.
		return randomStatistic <= testStatistic;
#endif
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
};

// Load set of genes from DB
QVector<Gene> loadGenes(QSqlDatabase &db, const QString &whereClause = "") {
	QVector<Gene> result;

	const QString sql = QString("SELECT l.Gene, x,y,z, SpeciesCount, Taxon FROM Loci l JOIN Conservation c ON l.Gene = c.Gene ORDER BY Chromosome, Start");
	QSqlQuery query(sql, db);
	while (query.next()) {
		Gene gene;
		gene.name = query.value(0).toString();
		gene.position.x = query.value(1).toDouble();
		gene.position.y = query.value(2).toDouble();
		gene.position.z = query.value(3).toDouble();
		gene.speciesCount = query.value(4).toInt();
		gene.taxon = query.value(5).toString();

		result.push_back(gene);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

#ifdef TAXON_TEST
using TaxonFrequencyMap = QMap<QString, double>;
TaxonFrequencyMap taxonFrequency;
#endif

// This implements the sphere test statistic for our particular test.
double sphereTestStatistic(const QVector<const Gene *> &genes) {
	if (genes.isEmpty())
		throw(QString("sphereTestStatistic(): Empty gene list provided!"));
	if (genes.size() == 1)
		return 0.0;


#ifdef TAXON_TEST
	// Compute frequencies in group
	TaxonFrequencyMap taxonFrequencyInGroup;
	for (const Gene *gene : genes) {
		taxonFrequencyInGroup[gene->taxon] += 1.0;
	}
	for (const QString &taxon : taxonFrequencyInGroup.keys()) {
		taxonFrequencyInGroup[taxon] /= (double)genes.size();
	}

	// For each taxon, calculate enrichment and keep the best
	double bestLogEnrichment = 0.0;
	for (const QString &taxon : taxonFrequencyInGroup.keys()) {
		const double enrichment = taxonFrequencyInGroup[taxon] / taxonFrequency[taxon];
		const double logEnrichment = log(enrichment);
		bestLogEnrichment = std::max(bestLogEnrichment, std::abs(logEnrichment));
	}

	return bestLogEnrichment;
#else
	// First, calculate average species count
	double averageSpeciesCount = 0.0;
	for (const Gene *gene : genes) {
		averageSpeciesCount += (double)gene->speciesCount;
	}
	averageSpeciesCount /= (double)genes.size();

	// Then variance
	double variance = 0.0;
	for (const Gene *gene : genes) {
		const double distanceFromAverage = (double)gene->speciesCount - averageSpeciesCount;
		variance += distanceFromAverage * distanceFromAverage;
	}

	return sqrt(variance);
#endif // TAXON_TEST
}

// This function performs the sphere test. An almost identical copy of this
// function exists in all similar sphere-test programs. This is a compromise
// between reusability and flexibility. Heavy parts of the procedure have been
// extracted to SphereTest.h.
void extractConservationFields(QSqlDatabase &db) {
	printf("Description of measured statistic:\n\t%s\n", statisticDescription);

	// Load genes
	const QVector<Gene> genes = loadGenes(db);
	printf("%d genes\n", genes.size());

#ifdef TAXON_TEST
	// Calculate general population taxon frequency
	for (const Gene &gene : genes) {
		taxonFrequency[gene.taxon] += 1.0;
	}
	printf("Base taxon frequencies:\n");
	for (const QString &taxon : taxonFrequency.keys()) {
		taxonFrequency[taxon] /= (double)genes.size();
		printf("\t%s: %.02f\n", taxon.toUtf8().data(), taxonFrequency[taxon]);
	}
#endif // TAXON_TEST

	QElapsedTimer timer;
	timer.start();

	// Generate a number of sphere samples

	printf("Generating %d sphere samples... ", sampleCount);
	int averageGenesInASphere = 0;
	QVector<WorkUnit<Gene>> workUnits =
		createWorkUnits(sphereRadius, minimumGeneCount, genes, sampleCount,
						&averageGenesInASphere);
	printf("Done.\n");
	printf("Average genes in a sphere: %d\n", averageGenesInASphere);

	// Then, for each of the samples, draw the same number of random samples of
	// the same gene count. Calculate the same metric and calculate a p-value.
	printf("Calculating p-values for %d sphere samples using %d random samples "
		   "for each... ",
		   sampleCount, sampleCount);
#pragma omp parallel for
	for (int i = 0; i < workUnits.size(); i++) {
		WorkUnit<Gene> &workUnit = workUnits[i];
		workUnit.calculatePValue(sampleCount);
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
	for (int i = 0; i < clusters.size(); i++) {
		printf("\tCluster %d: %d genes\n", i + 1, clusters[i].size());
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
		extractConservationFields(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
