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
sample random windows of 25 consecutive genes calculating Shannon entropy in the
window. The value of entropy will then be compared to the distribution of
entropy of random sets of 25 genes.
*/

#include <db/Communities.h>
#include <sampler/Entropy.h>
#include <utils/TsvReader.h>

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

void calculateEntropies(QSqlDatabase &db) {
	QVector<int> genes = db::loadGenes(db);

	printf("Pool of %d genes\n", genes.size());

	// Write entropies to file to study their distribution
	const std::string filename = "Results/CommunityEntropies.tsv";
	Tsv::Table table;
	Tsv::Row tableRow = {"SliceEntropy", "CloudEntropy"};
	table.push_back(tableRow);

	// Take continuous and random samples from the pool and calculate entropy.
	const int windowSize = 25;
	EntropySampler sampler(genes, windowSize);
	const int sampleCount = 100000;
	int chanceWinCount = 0;
	for (int i = 0; i < sampleCount; i++) {
		const double entropyOfSlice = sampler.sampleSlice();
		const double entropyOfCloud = sampler.sampleCloud();

		tableRow = {QString::number(entropyOfSlice).toStdString(),
					QString::number(entropyOfCloud).toStdString()};
		table.push_back(tableRow);

		if (entropyOfCloud <= entropyOfSlice)
			chanceWinCount++;
	}
	const double pValue = (double)chanceWinCount / (double)sampleCount;
	printf("p-value: %.03f\n", pValue);

	if (!Tsv::writeTsv(filename, table)) {
		throw(std::string("Failed to write output file: ") + filename);
	}
	printf("Written sampled entropies to file %s\n", filename.c_str());
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
