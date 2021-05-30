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
This program extracts the coexpression score distribution of "Continents" compartmentalization.
*/

#include "utils/PackedCoex.h"
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

struct Gene {
	QString name;
	Vec3D position;
	QString continent;

	// Index to square table of coexpression scores. This is both column and row index of the gene.
	// We use this to avoid searching gene pair coexpressions by name, thus improving performance.
	int coexIndex = 0;

	// One to rule them all
	static PackedCoex packedCoex;
};

PackedCoex Gene::packedCoex;

// Load set of genes from DB
QVector<Gene> loadGenes(QSqlDatabase &db, const QString &whereClause = "") {
	// Load packed coexpressions
	const QString coexFilename = QStringLiteral("Results/CoexPacked.bin");
	Gene::packedCoex.load(coexFilename);
	printf("Loaded packed coexpressions from file: %s\n",
		coexFilename.toUtf8().data());

	QVector<Gene> result;

	const QString sql = QString("SELECT l.Gene, x, y, z, Continent FROM Loci l JOIN Continents c ON l.Gene = c.Gene");
	QSqlQuery query(sql, db);
	while (query.next()) {
		Gene gene;
		gene.name = query.value(0).toString();
		gene.position.x = query.value(1).toDouble();
		gene.position.y = query.value(2).toDouble();
		gene.position.z = query.value(3).toDouble();
		gene.continent = query.value(4).toString();
		
		if (!Gene::packedCoex.geneToIndex.contains(gene.name)) {
			// Ignore genes where we don't know their coexpressions
			continue;
		}

		// Assign index to square table and we are done.
		gene.coexIndex = Gene::packedCoex.geneToIndex[gene.name];

		result.push_back(gene);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

void extractContinentStatistics(QSqlDatabase &db) {
	printf("Extracting coexpression score distribution of continents compartmentalization\n");

	// Load genes
	const QVector<Gene> genes = loadGenes(db);
	printf("%d genes\n", genes.size());

	// Split by continent
	QMap<QString, QVector<Gene>> continents;
	for (const Gene &gene : genes) {
		continents[gene.continent].push_back(gene);
	}
	for (const QString &continent : continents.keys()) {
		printf("\t%s:\t%d genes\n", continent.toUtf8().data(), continents[continent].size());
	}

	// Add a "Random" continent which is the entire genome
	continents.insert("Random", genes);

	const int sampleCount = 10000;

	// Write to file
	const QString &filename = "Results/ContinentCoexpressionScoreSamples.tsv";
	printf("Writing %d coexpression score samples to %s\n", sampleCount, filename.toUtf8().data());
	FILE *fp = fopen(filename.toUtf8().data(), "w");

	// Take a number of random samples for each continent.
	// One sample is one random pair of genes
	QStringList allContinentNames;
	for (const QString &continent : continents.keys()) {
		allContinentNames.push_back(continent);
	}
	fprintf(fp, "%s\n", allContinentNames.join("\t").toUtf8().data());
	
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, 10000);

	for (int i = 0; i < sampleCount; i++) {
		QStringList samples;
		for (const QString &continent : continents.keys()) {
			const QVector<Gene> &genes = continents[continent];
			const int geneIndexA = distribution(generator) % genes.size();
			int geneIndexB = geneIndexA;
			while (geneIndexB == geneIndexA) {
				geneIndexB = distribution(generator) % genes.size();
			}
			const Gene &geneA = genes[geneIndexA];
			const Gene &geneB = genes[geneIndexB];

			const double coexpressionScore = Gene::packedCoex.lookup(geneA.coexIndex, geneB.coexIndex) * 0.1;
			samples.push_back(QString::number(coexpressionScore));
		}

		fprintf(fp, "%s\n", samples.join("\t").toUtf8().data());
	} // end for (all samples)

	fclose(fp);

	// While we are here, let's also get some samples of pairs of genes, their proximity in space and their coexpression score.
	const QString proximityFilename = "Results/ProximityVsCoexpression.tsv";
	printf("Saving proximity and coexpression score samples to %s\n", proximityFilename.toUtf8().data());
	fp = fopen(proximityFilename.toUtf8().data(), "w");
	fprintf(fp, "Distance\tCoexpressionScore\n");
	for (int i = 0; i < sampleCount; i++) {
		const int geneIndexA = distribution(generator) % genes.size();
		int geneIndexB = geneIndexA;
		while (geneIndexB == geneIndexA) {
			geneIndexB = distribution(generator) % genes.size();
		}
		const Gene &geneA = genes[geneIndexA];
		const Gene &geneB = genes[geneIndexB];

		const double coexpressionScore = Gene::packedCoex.lookup(geneA.coexIndex, geneB.coexIndex) * 0.1;
		const double distance = Vec3D::distance(geneA.position, geneB.position);

		fprintf(fp, "%f\t%f\n", distance, coexpressionScore);
	} // end for (all samples)

	fclose(fp);
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
		extractContinentStatistics(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
