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
This module handles rendering chromosomes to SVG files. We use C++ templates so
that various types may be used as input. Chromosomes are represented by vectors
of "gene classes", where a gene class may be anything that can be used as a key
of a map. Strings and integers for instance. The elements are in order from
first to last gene of the chromosome and they classify the genes to a particular
class. The renderer will then output colored rectangles that show the position
of the classes within the chromosomes. We can output to both SVG for later
editing or HTML.
*/

#include <QFile>
#include <QMap>
#include <QString>
#include <QTextStream>
#include <QVector>

namespace Svg {

template <class ChromosomeName, class GeneClass>
void render(const QString &filename,
			const QMap<ChromosomeName, QVector<GeneClass>> &chromosomes,
			bool wrapToHtml) {

	using Chromosome = QVector<GeneClass>;

	// Count occurrences of each class so that we can auto-assign colors.
	QMap<GeneClass, int> classCounts;
	for (const Chromosome &chromosome : chromosomes.values()) {
		for (const GeneClass &geneClass : chromosome) {
			classCounts[geneClass]++;
		}
	}

	// Sort gene classes from least to most abundant
	QVector<GeneClass> geneClasses = classCounts.keys().toVector();
	std::sort(geneClasses.begin(), geneClasses.end(),
			  [&](const GeneClass &a, const GeneClass &b) {
				  return classCounts[a] < classCounts[b];
			  });

	// Assign colors. Black becomes the more abundant to improve contrast. Then
	// follows purple which is a mnemonic of "high frequency light", all the way
	// down to red which is low frequency. Then random colors.
	QVector<QString> colorPool = {"silver", "maroon",  "olive", "green",
								  "aqua",	"teal",	   "navy",	"fuchsia",
								  "red",	"lime",	   "blue",	"yellow",
								  "orange", "magenta", "black"};
	QMap<GeneClass, QString> classToColor;
	for (int i = geneClasses.size() - 1; i >= 0; i--) {
		if (colorPool.isEmpty()) {
			// Random color
			const int r = rand() % 256;
			const int g = rand() % 256;
			const int b = rand() % 256;
			const int rgb = (r << 16) | (g << 8) | b;
			QString color = QString::number(rgb, 16);
			while (color.size() < 6)
				color = QString("0") + color;
			classToColor[geneClasses[i]] = color;
		} else {
			classToColor[geneClasses[i]] = colorPool.back();
			colorPool.pop_back();
		}
	}

	// Space in pixels left empty between chromosomes
	const int padding = 25;
	const int chromosomeHeight = 20;

	// get maximum chromosome length
	int maxLength = 0;
	for (const Chromosome &c : chromosomes.values()) {
		maxLength = std::max(maxLength, (int)c.size());
	}

	printf("Saving to %s\n", filename.toUtf8().data());
	printf("%d chromosomes, max length: %d\n", chromosomes.count(), maxLength);

	const int imageWidth = 2 * padding + maxLength;
	int imageHeight =
		padding + (int)chromosomes.size() * (chromosomeHeight + padding);

	// Add height to fit a legend
	imageHeight += geneClasses.size() * padding;

	// Render to strings
	QVector<QString> svg;
	svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
	svg << QString("<svg width=\"%1\" height=\"%2\">")
			   .arg(imageWidth)
			   .arg(imageHeight);

	// Background gray box
	svg << QString("<rect x=\"%1\" y=\"%2\" width=\"%3\" "
				   "height=\"%4\" style=\"fill:%5\";\" />")
			   .arg(0)
			   .arg(0)
			   .arg(imageWidth)
			   .arg(imageHeight)
			   .arg("#777777");

	int x = padding;
	int y = padding;
	for (int i : chromosomes.keys()) {
		const Chromosome &chromosome = chromosomes[i];

		svg << "<g>";
		svg << QString("<text x=\"%1\" y=\"%2\">Chromosome %3</text>")
				   .arg(x)
				   .arg(y - 3)
				   .arg(i);

		for (GeneClass geneClass : chromosome) {
			svg << QString("<rect x=\"%1\" y=\"%2\" width=\"%3\" "
						   "height=\"%4\" style=\"fill:%5\";\" />")
					   .arg(x)
					   .arg(y)
					   .arg(1)
					   .arg(chromosomeHeight)
					   .arg(classToColor[geneClass]);

			x++;
		}
		svg << "</g>";

		x = padding;
		y += chromosomeHeight + padding;
	}

	// Legend
	for (GeneClass geneClass : geneClasses) {
		svg << QString("<text x=\"%1\" y=\"%2\" style=\"fill:%3\";\">%4</text>")
				   .arg(x)
				   .arg(y)
				   .arg(classToColor[geneClass])
				   .arg(geneClass);
		y += padding;
	}

	svg.push_back("</svg>");

	// Write to file
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		throw QString("Failed to open file %1 for writing\n").arg(filename);

	QTextStream out(&file);
	if (wrapToHtml)
		out << "<html><body><h1>Histone communities</h1>\n";
	for (const QString &line : svg) {
		out << line << "\n";
	}
	if (wrapToHtml)
		out << "</body></html>\n";
}

} // end namespace Svg