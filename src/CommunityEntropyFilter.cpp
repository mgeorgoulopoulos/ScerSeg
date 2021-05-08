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
This program will load genes and their "community" classification. Then it will
calculate Shannon entropy on a sliding window of 25 consecutive genes. Windows
will be filtered using a given threshold and then will be merged together where
they overlap. Finally, clusters of low entropy genes will be written to
database.
*/
#include <db/Communities.h>
#include <sampler/Entropy.h>
#include <utils/RenderSvg.h>

#include <QFile>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QString>
#include <QTextStream>
#include <QVector>

#include <algorithm>
#include <random>

namespace {

using Community = db::Community;

// Combines the overlapping windows and reports the result in console.
void reportTightGroups(const QVector<bool> &tightGroups,
					   const QVector<int> &genes) {
	using StartEnd = QPair<int, int>;
	QVector<StartEnd> groups;

	if (tightGroups.size() != genes.size())
		throw(QString("Number of thresholded entropies (%1) does not match "
					  "number of genes (%2)")
				  .arg(tightGroups.size())
				  .arg(genes.size()));

	int counter = 0;

	StartEnd currentGroup;
	bool halfComplete = false;
	for (int i = 1; i < tightGroups.size(); i++) {
		const bool a = tightGroups[i - 1];
		const bool b = tightGroups[i];

		if (a)
			counter++;

		if (a == b)
			continue;

		if (b) {
			// Start of a region
			currentGroup.first = i;
			halfComplete = true;
			continue;
		}

		// End of a region
		currentGroup.second = i;
		groups.push_back(currentGroup);
		halfComplete = false;
	}

	printf("COUNTER: %d\n", counter);

	// Handle edge case
	if (halfComplete) {
		currentGroup.second = tightGroups.size() - 1;
		groups.push_back(currentGroup);
		halfComplete = false;
	}

	int geneCount = 0;
	for (const StartEnd &group : groups) {
		printf("%d\t%d (%d genes): ", group.first, group.second,
			   group.second - group.first + 1);

		// Also tell which clusters are contained
		for (int i = group.first; i <= group.second; i++) {
			printf("%d ", genes[i]);
		}
		printf("\n");

		geneCount += group.second - group.first + 1;
	}
	printf("\n%d genes total.\n", geneCount);
}

void execNonQuery(QSqlDatabase &db, const QString &sql) {
	QSqlQuery query = db.exec(sql);
	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));
}

struct Gene {
	QString name;
	int chromosome = 0;
	int community = 0;
};

QVector<Gene> loadGenes(QSqlDatabase &db) {
	QVector<Gene> result;

	const QString sql = "SELECT l.Gene, Chromosome, Community FROM Communities "
						"c LEFT JOIN Loci l ON "
						"c.Gene = l.Gene ORDER BY Chromosome, Start";
	QSqlQuery query(sql, db);

	while (query.next()) {
		Gene gene;
		gene.name = query.value(0).toString();
		gene.chromosome = query.value(1).toInt();
		gene.community = query.value(2).toInt();

		result.push_back(gene);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

// Writes the filtered (<bool>=true) genes to database
void populateTightCommunities(QSqlDatabase &db,
							  const QVector<bool> &tightGroups) {
	// Read gene names in the same order they appear in the state list.
	QVector<Gene> genes = loadGenes(db);

	if (genes.size() != tightGroups.size())
		throw(QString("Gene names (%d) and tight groups (%d): size mismatch")
				  .arg(genes.size())
				  .arg(tightGroups.size()));

	QVector<QString> tightGenes;
	for (int i = 0; i < genes.size(); i++) {
		if (!tightGroups[i])
			continue;

		tightGenes.push_back(genes[i].name);
	}

	// Recreate and populate TightComminities table.
	execNonQuery(db, "DROP TABLE IF EXISTS TightCommunities");
	execNonQuery(db, "CREATE TABLE TightCommunities(Gene TEXT PRIMARY KEY)");

	const QString sqlInsert =
		"INSERT INTO TightCommunities (Gene) VALUES (:Gene)";
	QSqlQuery queryInsert(db);
	if (!queryInsert.prepare(sqlInsert))
		throw QString("Failed to create query: %1").arg(sqlInsert);

	for (const QString &gene : tightGenes) {
		queryInsert.bindValue(":Gene", gene);
		if (!queryInsert.exec())
			throw QString("Failed to exec query: %1").arg(sqlInsert);
	}
}

// Renders to SVG
void renderSvg(QSqlDatabase &db, const QVector<bool> &tightGroups) {

	QVector<Gene> genes = loadGenes(db);
	if (genes.size() != tightGroups.size())
		throw(QString("Gene names (%d) and tight groups (%d): size mismatch")
				  .arg(genes.size())
				  .arg(tightGroups.size()));

	// Convert to data structure expected by svg renderer
	using ChromosomeName = int;
	using GeneClass = int;
	QMap<ChromosomeName, QVector<GeneClass>> renderData;

	for (int i = 0; i < genes.size(); i++) {
		const GeneClass geneClass = tightGroups[i] ? genes[i].community : 0;
		renderData[genes[i].chromosome].push_back(geneClass);
	}

	// Render to both SVG and HTML
	Svg::render("Results/TightCommunities.svg", renderData, false);
	Svg::render("Results/TightCommunities.html", renderData, true);
}

void calculateEntropies(QSqlDatabase &db) {
	QVector<int> genes = db::loadGenes(db);

	printf("Pool of %d genes\n", genes.size());

	// Entropy threshold. This value correspond of a (non-adjusted) p-value of
	// 0.001
	const double entropyThreshold = 1.28;
	QVector<bool> tightGroups;
	tightGroups.resize(genes.size());
	tightGroups.fill(false);

	// Measure entropy on a sliding window
	const int windowSize = 25;
	EntropySampler sampler(genes, windowSize);
	int valuesUnderThreshold = 0;

	// Do the sliding window thing
	const int sampleCount = genes.size() - windowSize;
	for (int i = 0; i < sampleCount; i++) {
		// Calculate entropy of slice
		const double sliceEntropy = sampler.sampleSlice(i);

		// Mark tight groups
		if (sliceEntropy < entropyThreshold) {
			valuesUnderThreshold++;

			// mark the entire window as "tight"
			for (int t = 0; t < windowSize; t++) {
				tightGroups[i + t] = true;
			}
		}
	}

	printf("Values under threhsold (%.02f): %d\n\n", entropyThreshold,
		   valuesUnderThreshold);

	// Maybe some tight groups overlap - combine to (start,end) pairs and
	// report.
	reportTightGroups(tightGroups, genes);

	// Export genes for later use
	populateTightCommunities(db, tightGroups);

	// Render linear chromosomes
	renderSvg(db, tightGroups);
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
		calculateEntropies(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
