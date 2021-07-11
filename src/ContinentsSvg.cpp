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
This program reads "continent" classification from database and renders an SVG
image with the chromosomes and areas of them classified to color-coded
continents.
*/

#include "Utils/RenderSvg.h"

#include <QFile>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QString>
#include <QTextStream>
#include <QVector>

#include <algorithm>

namespace {
using Continent = QString;
using Chromosome = QVector<Continent>;

QMap<int, Chromosome> loadChromosomes(QSqlDatabase &db) {
	QMap<int, Chromosome> result;

	const QString sql =
		"SELECT Chromosome, Field FROM Loci l JOIN ContinentFields c ON l.Gene = c.Gene ORDER BY Chromosome, Start";
	QSqlQuery query(sql, db);
	while (query.next()) {
		const int chromosome = query.value(0).toInt();
		const Continent continent = query.value(1).toString();
		result[chromosome].push_back(continent);
	}

	if (query.lastError().type() != QSqlError::NoError)
		throw(QString("Failed to process query: %1\nDBTEXT: %2")
				  .arg(sql)
				  .arg(query.lastError().databaseText()));

	return result;
}

void renderContinents(QSqlDatabase &db) {
	QMap<int, Chromosome> chromosomes = loadChromosomes(db);

	// Assign colors
	QMap<Continent, QString> continentColor;
	continentColor.insert("Tethys", "black");
	continentColor.insert("Laurasia", "orange");
	continentColor.insert("Godwana", "fuchsia");
	continentColor.insert("Antarctica", "red");

	Svg::render("Results/Continents.html", chromosomes, true, continentColor);
	Svg::render("Results/Continents.svg", chromosomes, false, continentColor);
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
		renderContinents(db);
	} catch (QString errorMessage) {
		printf("ERROR: %s\n", errorMessage.toUtf8().data());
		return 0;
	}

	db.close();

	printf("Full success\n");

	return 0;
}
