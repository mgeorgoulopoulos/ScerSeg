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
This program samples distances in 3D space between: (1) genes found in
TightCommunities db table and (2) random pairs of genes. Also, we always sample pairs
from different chromosomes in order to remove the community advantage of being
consecutive. Results are written to TSV file (one column for (1) and a second
for (2)) for further analysis.
*/

#include <utils/Vec3D.h>

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

struct Gene {
	QString name;
	int chromosome = 0;
	Vec3D position;
};

QVector<Gene> loadGenes(QSqlDatabase &db, const QString &whereClause = "") {
	QVector<Gene> result;

	const QString sql = QString("SELECT Gene, x, y, z, Chromosome FROM Loci %1 "
								"ORDER BY Chromosome, Start")
							.arg(whereClause);
	QSqlQuery query(sql, db);
	while (query.next()) {
		Gene gene;
		gene.name = query.value(0).toString();
		gene.position.x = query.value(1).toDouble();
		gene.position.y = query.value(2).toDouble();
		gene.position.z = query.value(3).toDouble();
		gene.chromosome = query.value(4).toInt();

		result.push_back(gene);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

void sampleDistances(QSqlDatabase &db) {
	const QVector<Gene> genes = loadGenes(db);
	const QVector<Gene> tightGenes =
		loadGenes(db, "WHERE Gene IN (SELECT Gene FROM TightCommunities)");

	printf("%d genes\n", genes.size());
	printf("%d genes in tight groups\n", tightGenes.size());

	// Break into chromosomes
	QMap<int, QVector<Gene>> tightGroupsMap;
	for (const Gene &gene : tightGenes) {
		tightGroupsMap[gene.chromosome].push_back(gene);
	}

	QVector<QVector<Gene>> tightGroups;
	for (const int chromosome : tightGroupsMap.keys()) {
		printf("Chromosome %d: group of size %d\n", chromosome,
			   tightGroupsMap[chromosome].size());
		tightGroups.push_back(tightGroupsMap[chromosome]);
	}

	// Do the same for gene pool
	QMap<int, QVector<Gene>> geneMap;
	for (const Gene &gene : genes) {
		geneMap[gene.chromosome].push_back(gene);
	}

	QVector<QVector<Gene>> geneGroups;
	for (const int chromosome : geneMap.keys()) {
		geneGroups.push_back(geneMap[chromosome]);
	}

	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, genes.size() - 1);
	auto randomInt = [&](int maximumExcl) {
		return distribution(generator) % maximumExcl;
	};

	// Output distributions to tsv
	const QString filename = "Results/TightCommunityDistances.tsv";
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		throw QString("Failed to open file %1 for writing\n").arg(filename);
	QTextStream out(&file);

	// Take samples
	int chanceWinCount = 0;
	const int sampleCount = 1000;
	int atomicSampleCount = 0;
	for (int i = 0; i < sampleCount; i++) {
		// Visit all pairs of groups
		for (int groupIndex1 = 0; groupIndex1 < tightGroups.size();
			 groupIndex1++) {
			const QVector<Gene> &group1 = tightGroups[groupIndex1];
			const QVector<Gene> &randomGroup1 = geneGroups[groupIndex1];
			for (int groupIndex2 = groupIndex1 + 1;
				 groupIndex2 < tightGroups.size(); groupIndex2++) {
				const QVector<Gene> &group2 = tightGroups[groupIndex2];
				const QVector<Gene> &randomGroup2= geneGroups[groupIndex2];

				// Take random genes from group1 & 2
				const Gene &tightGene1 = group1[randomInt(group1.size())];
				const Gene &tightGene2 = group2[randomInt(group2.size())];
				const double tightDistance =
					Vec3D::distance(tightGene1.position, tightGene2.position);

				// Take random genes from general population
				const Gene &randomGene1 = randomGroup1[randomInt(randomGroup1.size())];
				const Gene &randomGene2 = randomGroup2[randomInt(randomGroup2.size())];
				const double randomDistance =
					Vec3D::distance(randomGene1.position, randomGene2.position);

				out << tightDistance << '\t' << randomDistance << '\n';

				if (tightDistance < randomDistance)
					chanceWinCount++;
				atomicSampleCount++;
			}
		}
	} // end for (N samples)

	double pValue = (double)chanceWinCount / (double)atomicSampleCount;

	printf("%d samples taken\n", atomicSampleCount);
	printf("p-value: %f\n", pValue);
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
		sampleDistances(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
