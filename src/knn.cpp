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
This program reads a scaffold of genes from database and then runs k-nearest
neighbors to mount all other genes to the scaffold. Distances are measured in
histone space. Result of the classification is then saved to Communities
database table.
*/

#define HISTONE_COUNT 9

#include <QFile>
#include <QMap>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QString>
#include <QVariant>
#include <QVector>

#include <algorithm>

namespace {
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

QVector<Gene> loadGenes(QSqlDatabase &db) {
	QVector<Gene> result;

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

		result.push_back(gene);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

using GeneToCommunity = QMap<QString, int>;
GeneToCommunity loadScaffold(QSqlDatabase &db) {
	GeneToCommunity result;

	const QString sql = "SELECT Gene, Community FROM CommunityCenters";
	QSqlQuery query(sql, db);

	while (query.next()) {
		const QString gene = query.value(0).toString();
		const int community = query.value(1).toInt();
		result[gene] = community;
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

// Classifies all genes based on scaffold
GeneToCommunity knnClassify(const QVector<Gene> &genes,
							const GeneToCommunity &scaffold) {
	const int k = 3;

	printf("Classifying %d genes using a scaffold of %d genes\n", genes.size(),
		   scaffold.size());
	printf("Using:\n\tk=%d\n", k);

	// Our result will be a superset of the scaffold
	GeneToCommunity result = scaffold;

	// Get the subset of scaffold genes
	QVector<Gene> scaffoldGenes;
	scaffoldGenes.reserve(scaffold.size());
	for (const Gene &gene : genes) {
		if (scaffold.contains(gene.name)) {
			scaffoldGenes.push_back(gene);
		}
	}

	if (scaffoldGenes.size() != scaffold.size())
		throw(QString("Something went wrong in getting the scaffold gene "
					  "subset. Got %1 genes instead of %2")
				  .arg(scaffoldGenes.size())
				  .arg(scaffold.size()));

	// find k nearest neighbours
	QVector<Gene> knn;
	for (const Gene &gene : genes) {
		auto closest = [&](const Gene &a, const Gene &b) {
			return gene.distance(a) < gene.distance(b);
		};

		// Ignore genes already in the scaffold
		if (scaffold.contains(gene.name))
			continue;

		// Find the k nearest neighbors
		knn.resize(0);
		for (const Gene &scaffoldGene : scaffoldGenes) {
			knn.push_back(scaffoldGene);
			std::sort(knn.begin(), knn.end(), closest);
			if (knn.size() > k)
				knn.resize(k);
		}

		// Count votes (occurrences of each community in the k neighbors)
		using CommunityToVoteCount = QMap<int, int>;
		CommunityToVoteCount votes;
		int maxVotedCommunity = -1;
		int maxVotes = 0;
		for (const Gene &g : knn) {
			const int community = scaffold[g.name];
			int &v = votes[community];
			v++;
			if (maxVotedCommunity < 0 || v > maxVotes) {
				maxVotedCommunity = community;
				maxVotes = v;
			}
		}
		if (maxVotedCommunity < 0) {
			throw(QString("Could not classify gene: %1").arg(gene.name));
		}

		// Save result
		result[gene.name] = maxVotedCommunity;
	}

	return result;
}

void writeCommunities(QSqlDatabase &db, const GeneToCommunity &communities) {
	const QString sqlDrop = "DROP TABLE IF EXISTS Communities";
	QSqlQuery queryDrop(db);
	if (!queryDrop.prepare(sqlDrop))
		throw QString("Failed to crete query: %1").arg(sqlDrop);
	if (!queryDrop.exec())
		throw QString("Failed to exec query: %1").arg(sqlDrop);

	const QString sqlCreate =
		"CREATE TABLE Communities(Gene TEXT PRIMARY KEY, Community INTEGER)";
	QSqlQuery queryCreate(db);
	if (!queryCreate.prepare(sqlCreate))
		throw QString("Failed to crete query: %1").arg(sqlCreate);
	if (!queryCreate.exec())
		throw QString("Failed to exec query: %1").arg(sqlCreate);

	// Insert data
	const QString sqlInsert =
		"INSERT INTO Communities(Gene, Community) VALUES (:Gene, :Community)";
	QSqlQuery queryInsert(db);
	if (!queryInsert.prepare(sqlInsert))
		throw QString("Failed to create query: %1").arg(sqlInsert);

	for (const QString &gene : communities.keys()) {
		queryInsert.bindValue(":Gene", gene);
		queryInsert.bindValue(":Community", communities[gene]);
		if (!queryInsert.exec())
			throw QString("Failed to exec query: %1").arg(sqlInsert);
	}

	printf("Created table 'Communities' with %d rows\n", communities.size());
}

void classifyAndUpdateDatabase(QSqlDatabase &db) {
	// Load genes
	QVector<Gene> genes = loadGenes(db);

	// Load scaffold
	GeneToCommunity scaffold = loadScaffold(db);

	// Apply KNN classification
	GeneToCommunity communities = knnClassify(genes, scaffold);

	// Update db
	writeCommunities(db, communities);
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
		classifyAndUpdateDatabase(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
