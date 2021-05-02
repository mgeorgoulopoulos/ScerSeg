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
This program loads histone modifications from PromoterHistonesPatched table and
then, for each gene: it calculates the gene's distance in histone space against
all other genes (as average Euclidean distance) and the same average distance in
its local group of 5(the gene plus 2 genes from both its left and right).Then,
the "community score" is calculated as averageDistanceFromAll /
averageDistanceInLocalGroup. The scores are then inserted to "CommunityScores"
database table.
*/

#include <QMap>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QVariant>
#include <QVector>

#define HISTONE_COUNT 9

struct Gene {
	QString name;
	int chromosome = 1;
	double histones[HISTONE_COUNT];
	double score = 0.0;

	QString tsvString() const {
		QString result = QString("%1\t%2").arg(name).arg(chromosome);
		for (int i = 0; i < HISTONE_COUNT; i++) {
			result += QString("\t%1").arg(histones[i]);
		}
		return result;
	}

	double distance(const Gene &other) const {
		double accum = 0.0;
		for (int i = 0; i < HISTONE_COUNT; i++) {
			const double d = other.histones[i] - histones[i];
			accum += d * d;
		}
		return sqrt(accum);
	}
};

// Calculates average histone distance in a gene group (everything vs
// everything).
double averageDistance(const QVector<Gene> &group) {
	if (group.size() <= 1)
		return 0.0;
	double averageDistance = 0.0;
	int count = 0;
	for (int i = 0; i < group.size(); i++) {
		for (int j = i + 1; j < group.size(); j++) {
			averageDistance += group[i].distance(group[j]);
			count++;
		}
	}

	Q_ASSERT(count > 0);
	averageDistance /= (double)count;

	return averageDistance;
}

// Average distances are calculated per-chromosome. We visit all genes and
// calculate each one's score. We save the score in the gene struct itself.
void processChromosome(int chromosomeNumber, QVector<Gene> &genes) {
	printf("Processing: chromosome %d: %d genes\n", chromosomeNumber,
		   genes.size());
	if (genes.isEmpty()) {
		// Nothing to do here
		return;
	}

	// Sliding window of 5 genes
	QVector<Gene> localGroup;
	for (int i = 0; i < genes.size(); i++) {
		Gene &gene = genes[i];

		// Find distance of gene from all other genes
		double globalAverageDitsance = 0.0;
		for (int j = 0; j < genes.size(); j++) {
			globalAverageDitsance += gene.distance(genes[j]);
		}
		globalAverageDitsance /= (double)genes.size();

		// Then, find average distance of the local group of 5 genes
		localGroup.resize(0);
		for (int j = i - 2; j <= i + 2; j++) {
			if (j < 0 || j >= genes.size())
				continue;
			localGroup.push_back(genes[j]);
		}

		double localGroupAverageDistance = averageDistance(localGroup);

		// Prevent division by zero by clamping to a low value. This value is
		// average minus 3 standard deviations.
		const double veryLowValue = 0.17;
		localGroupAverageDistance =
			std::max(veryLowValue, localGroupAverageDistance);

		// Save the score to gene
		gene.score = globalAverageDitsance / localGroupAverageDistance;
	}
}

// Deletes ard re-creats CommunityScores table, which assigns a score to each
// key (Gene name)
void writeScores(QSqlDatabase &db,
				 const QMap<int, QVector<Gene>> &chromosomes) {
	const QString sqlDrop = "DROP TABLE IF EXISTS CommunityScores";
	QSqlQuery queryDrop(db);
	if (!queryDrop.prepare(sqlDrop))
		throw QString("Failed to crete query: %1").arg(sqlDrop);
	if (!queryDrop.exec())
		throw QString("Failed to exec query: %1").arg(sqlDrop);

	const QString sqlCreate =
		"CREATE TABLE CommunityScores(Gene TEXT PRIMARY KEY, Score REAL)";
	QSqlQuery queryCreate(db);
	if (!queryCreate.prepare(sqlCreate))
		throw QString("Failed to crete query: %1").arg(sqlCreate);
	if (!queryCreate.exec())
		throw QString("Failed to exec query: %1").arg(sqlCreate);

	// Insert data
	const QString sqlInsert =
		"INSERT INTO CommunityScores (Gene, Score) VALUES (:Gene, :Score)";
	QSqlQuery queryInsert(db);
	if (!queryInsert.prepare(sqlInsert))
		throw QString("Failed to create query: %1").arg(sqlInsert);

	for (const int c : chromosomes.keys()) {
		const QVector<Gene> &genes = chromosomes[c];
		for (const Gene &gene : genes) {
			queryInsert.bindValue(":Gene", gene.name);
			queryInsert.bindValue(":Score", gene.score);
			if (!queryInsert.exec())
				throw QString("Failed to exec query: %1").arg(sqlInsert);
		} // end for (all genes)
	}	  // end for (all chromosomes)
}

// Loads genes from db
QMap<int, QVector<Gene>> loadGenes(QSqlDatabase &db) {
	QMap<int, QVector<Gene>> result;

	const QString sql =
		"SELECT Chromosome, h.* FROM HistonesPromoterPatched h JOIN Loci l ON "
		"h.Gene = l.Gene ORDER BY Chromosome, Start";
	QSqlQuery query(sql, db);
	while (query.next()) {
		Gene gene;
		gene.chromosome = query.value(0).toInt();
		gene.name = query.value(1).toString();
		for (int i = 0; i < HISTONE_COUNT; i++) {
			gene.histones[i] = query.value(i + 2).toDouble();
		}

		result[gene.chromosome].push_back(gene);
	}

	return result;
}

void calculateScoresAndUpdateDatabase(QSqlDatabase &db) {
	// Load genes
	QMap<int, QVector<Gene>> chromosomes = loadGenes(db);

	// process all chromosomes
	for (const int c : chromosomes.keys()) {
		QVector<Gene> &genes = chromosomes[c];
		processChromosome(c, genes);
	}

	// Update db
	writeScores(db, chromosomes);
}

int main(int argc, char *argv[]) {

	const QString &filename = "Results/yeast.sqlite";
	QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE");
	db.setDatabaseName(filename);
	if (!db.open()) {
		printf("Failed to open file: %s\n", filename.toUtf8().data());
	}

	try {
		calculateScoresAndUpdateDatabase(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
