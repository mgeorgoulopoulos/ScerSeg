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
This module defines the PackedCoex class, which is responsible to serialize /
deserialize the binary coexpression format. We use this format for speed during
the respective sphere test.
*/

#include <QMap>
#include <QString>
#include <QVector>

struct PackedCoex {
	QVector<QString> genes;
	QMap<QString, int> geneToIndex;
	QVector<unsigned char> coex;

	void load(const QString &filename) {
		genes.clear();
		coex.clear();

		FILE *fp = fopen(filename.toUtf8().data(), "rb");
		QString gene;
		while (true) {
			if (feof(fp)) {
				throw(QString("Packed file is corrupted. Haven't finished "
							  "loading gene names. I have only %1")
						  .arg(genes.size()));
			}

			const char c = getc(fp);
			if (c == 0) {
				// null-terminated
				if (gene.isEmpty()) {
					break;
				} else {
					genes.push_back(gene);
					gene.clear();
				}
				continue;
			} // end if null

			gene += c;
		}

		// Create index
		for (int i = 0; i < genes.size(); i++) {
			geneToIndex.insert(genes[i], i);
		}

		// Then load coexpression scores
		coex.reserve(genes.size() * genes.size());
		for (int i = 0; i < genes.size(); i++) {
			for (int j = 0; j < genes.size(); j++) {
				if (feof(fp)) {
					throw(QString("Packed file is too short. I need %1 values "
								  "but I have %d")
							  .arg(genes.size() * genes.size())
							  .arg(coex.size()));
				}
				coex.push_back((unsigned char)getc(fp));
			}
		}

		fclose(fp);
	}

	unsigned char lookup(int gene1, int gene2) {
		return coex[gene1 * genes.size() + gene2];
	}

	unsigned char lookup(const QString &gene1, const QString &gene2) {
		int i = geneToIndex[gene1];
		int j = geneToIndex[gene2];
		return lookup(i, j);
	}

	void validate(const QString &gene1, const QString &gene2,
				  const unsigned char expectedScore) {
		const unsigned char lookupScore = lookup(gene1, gene2);
		if (expectedScore != lookupScore)
			throw(QString("Validation failed for %1 -> %2. Expected %3, got %4")
					  .arg(gene1)
					  .arg(gene2)
					  .arg(expectedScore)
					  .arg(lookupScore));
	}
};