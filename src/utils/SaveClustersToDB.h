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
This module provides a function to save a number of clusters to an sqlite
database. We use a naming scheme which assigns the letters of the English
alphabet in order of increasing size. So 'A' is the smallest and so on. Our
convention is to name these clusters "Fields". In a field we expect one or a few
species of plants to grow. The clusters resulting from this study are selected
in such a way that they exhibit a statistically significant preference of one
thing per cluster. A certain histone modification profile, for instance. Thus,
the newly created sqlite table has a "Gene" column and a "Field" column
containing the assignment of the gene to a particular field ('A', 'B' and such).
provided.
*/

#ifndef _SAVE_CLUSTERS_TO_DB_H_
#define _SAVE_CLUSTERS_TO_DB_H_

#include <QSet>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QString>
#include <QVariant>
#include <QVector>

namespace db {

void writeClusters(const QVector<QSet<QString>> &geneClusters, QSqlDatabase &db,
				   const QString &tableName) {

	// Drop possibly pre-existing instance of the table
	const QString sqlDrop = QString("DROP TABLE IF EXISTS %1").arg(tableName);
	QSqlQuery queryDrop(db);
	if (!queryDrop.prepare(sqlDrop))
		throw QString("Failed to crete query: %1").arg(sqlDrop);
	if (!queryDrop.exec())
		throw QString("Failed to exec query: %1").arg(sqlDrop);

	// (Re)-create the table
	const QString sqlCreate =
		QString("CREATE TABLE %1 (Gene TEXT PRIMARY KEY, Field TEXT)")
			.arg(tableName);
	QSqlQuery queryCreate(db);
	if (!queryCreate.prepare(sqlCreate))
		throw QString("Failed to crete query: %1").arg(sqlCreate);
	if (!queryCreate.exec())
		throw QString("Failed to exec query: %1").arg(sqlCreate);

	// Insert data
	const QString sqlInsert =
		QString("INSERT INTO %1(Gene, Field) VALUES (:Gene, :Field)")
			.arg(tableName);
	QSqlQuery queryInsert(db);
	if (!queryInsert.prepare(sqlInsert))
		throw QString("Failed to create query: %1").arg(sqlInsert);

	for (int i = 0; i < geneClusters.size(); i++) {
		QString fieldName = 'A' + i; // Name the fields as 'A', 'B', 'C' ...
		const QSet<QString> &genes = geneClusters[i];
		for (const QString &gene : genes) {
			queryInsert.bindValue(":Gene", gene);
			queryInsert.bindValue(":Field", fieldName);
			if (!queryInsert.exec())
				throw QString("Failed to exec query: %1").arg(sqlInsert);
		} // end for (all genes)
	}	  // end for (all clusters)
}

} // end namespace db

#endif // _SAVE_CLUSTERS_TO_DB_H_
