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
This module handles database access to Communities table. This table assigns a
numerical class to each gene, placing it in one of the histone classes.
*/

#include <QSqlDatabase>
#include <QSqlQuery>
#include <QString>
#include <QVariant>
#include <QVector>
#include <QSqlError>

namespace db {

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

} // end namespace db
