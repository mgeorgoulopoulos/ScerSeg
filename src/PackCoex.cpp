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
This program preprocesses the massive zip file with the SGDgene files. One for 
each gene, and extracts coexpression scores in a square byte matrix laid out 
in one single array of bytes. This reduces the multiple gigabytes of data to 
~35 MB, which allows us then to perform the heavy sphere-test on the data.
The conversion to bytes is lossless. Floating point scores provided range 
from 0.0 to 25.0, always specifying a single decimal digit.
We also provide rudimentary testing, against 20 random hand-picked records.
*/

#include "utils/PackedCoex.h"

#include <QApplication>

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

struct CoexRecord {
	QString gene1;
	QString gene2;
	double score = 0.0;
};

QVector<CoexRecord> load(QSqlDatabase &db) {
	QVector<CoexRecord> result;

	const QString sql = QString("SELECT Gene1, Gene2, Score FROM Coex");
	QSqlQuery query(sql, db);

	QSet<QString> genes;

	while (query.next()) {
		CoexRecord record;
		record.gene1 = query.value(0).toString();
		record.gene2 = query.value(1).toString();
		record.score = query.value(2).toDouble();

		result.push_back(record);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

#define CREATE_PACK

void packCoex(QSqlDatabase &db) {
	const QString packedFilename = QStringLiteral("CoexPacked.bin");
#ifdef CREATE_PACK
	printf("Loading coexpression records from DB... ");
	QVector<CoexRecord> records = load(db);
	printf("Done\n%d records\n", records.size());

	// Convert to map for easy access.
	printf("Constructing lookup matrix... ");
	QSet<QString> allGenes;
	QMap<QString, QMap<QString, double>> matrix;
	for (const CoexRecord &record : records) {
		matrix[record.gene1][record.gene2] = record.score;
		allGenes.insert(record.gene1);
		allGenes.insert(record.gene2);
	}
	printf("Done\n");

	// Create index
	QVector<QString> genes = allGenes.toList().toVector();
	std::sort(genes.begin(), genes.end());

	// Write to binary file
	printf("Writing binary file %s ... ", packedFilename.toUtf8().data());
	FILE *fp = fopen(packedFilename.toUtf8().data(), "wb");
	for (const QString &gene : genes) {
		fprintf(fp, "%s", gene.toUtf8().data());
		putc(0, fp); // null-terminate
	}
	putc(0, fp); // null-terminate all strings
	for (const QString &gene1 : genes) {
		for (const QString &gene2 : genes) {
			const double score = matrix[gene1][gene2];

			// pack to 1 byte per score - zero data loss as scores are up
			// to 21.0 with only one decimal digit
			const unsigned char packed = (unsigned char)(score * 10);
			putc(packed, fp);
		}
	}
	fclose(fp);
	printf("Done\n");
#endif
	// Validate exported binary file
	printf("Validating packed file...\n");
	PackedCoex packed;
	packed.load(packedFilename);
	printf("\tLoading successful\n");

#ifdef CREATE_PACK
	// Verify gene count
	QSet<QString> loadedGeneSet = packed.genes.toList().toSet();
	if (matrix.keys().toSet() != loadedGeneSet) {
		throw(QString("Different set of genes loaded from packed file"));
	}
	printf("\tGene sets identical\n");
#endif

	// Verify hand-picked records
	packed.validate("YCR107W", "YJL101C", 17);
	packed.validate("YCR107W", "YPL280W", 12);
	packed.validate("YCR107W", "YDR330W", 11);
	packed.validate("YCR107W", "YMR301C", 11);
	packed.validate("YKR074W", "YPR020W", 13);
	packed.validate("YKR074W", "YFR024C", 2);
	packed.validate("YKR074W", "YNR042W", 5);
	packed.validate("YKR074W", "YOR391C", 5);
	packed.validate("YDR118W", "YHR216W", 9);
	packed.validate("YDR118W", "YDR368W", 11);
	packed.validate("YDR118W", "YGL251C", 13);
	packed.validate("YJL081C", "YGL169W", 15);
	packed.validate("YJL081C", "YLR129W", 14);
	packed.validate("YJL081C", "YDR432W", 14);
	packed.validate("YJL081C", "YNL309W", 11);
	packed.validate("YBL078C", "YHL023C", 12);
	packed.validate("YBL078C", "YLR339C", 2);
	packed.validate("YBL078C", "YML053C", 9);
	packed.validate("YBL078C", "YHR068W", 6);

	// This one we found the hard way:
	packed.validate("YMR312W", "YHL011C", 27);

	printf("\tHand-picked validation successful\n");
}

} // end anonymous namespace

int main(int argc, char *argv[]) {
	const QString &filename = "Results/coex.sqlite";

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
		packCoex(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
