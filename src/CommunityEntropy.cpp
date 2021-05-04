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

using Community = int;

QVector<Community> loadGenes(QSqlDatabase &db) {
	QVector<Community> result;

	const QString sql = "SELECT Community FROM Communities c JOIN Loci l ON "
						"c.Gene = l.Gene ORDER BY Chromosome, Start";
	QSqlQuery query(sql, db);
	while (query.next()) {
		const Community community = query.value(0).toInt();

		result.push_back(community);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

// Helper class for sampling (uniformly) random genes and calculating Shannon
// entropy
class EntropySampler {
  public:
	EntropySampler(const QVector<Community> &pool, int windowSize = 25)
		: pool(pool), cloudDistribution(0, pool.size() - 1),
		  poolDistribution(0, pool.size() - windowSize) {
		window.resize(windowSize);
		if (!pool.empty()) {
			int maxValue = *std::max_element(pool.begin(), pool.end());
			counts.resize(maxValue + 1);
			counts.fill(0);
		}
	}

	// Entropy of a random continuous slice of 'windowSize'
	double sampleSlice() const {
		const int start = poolDistribution(generator);
		for (int i = 0; i < window.size(); i++) {
			window[i] = pool[start + i];
		}
		return entropy(window);
	}

	// Entropy of a random set of genes from the pool
	double sampleCloud() const {
		for (int i = 0; i < window.size(); i++) {
			window[i] = pool[cloudDistribution(generator)];
		}
		return entropy(window);
	}

	double entropy(const QVector<Community> &sample) const {
		counts.fill(0);
		for (const int v : sample) {
			Q_ASSERT(v >= 0 && v < counts.size());
			counts[v]++;
		}

		// calculate Shannon entropy
		double e = 0.0;
		for (const int c : counts) {
			if (c <= 0) {
				// zero probability
				continue;
			}
			// probability of the symbol
			const double p = (double)c / (double)sample.size();

			// Shannon formula
			e -= p * log2(p);
		}

		return e;
	}

  private:
	mutable QVector<int> counts;
	mutable QVector<Community> window;
	const QVector<Community> &pool;
	mutable std::default_random_engine generator;
	mutable std::uniform_int_distribution<int> cloudDistribution;
	mutable std::uniform_int_distribution<int> poolDistribution;
};

void calculateEntropies(QSqlDatabase &db) {
	QVector<int> genes = loadGenes(db);

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
