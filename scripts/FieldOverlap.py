# Copyright 2021 Michael Georgoulopoulos mgeorgoulopoulos@gmail.com

# This code is licensed under the MIT license. See LICENSE for details.

# Calculate overlap of all extracted clusters.


import sqlite3

tables = ["PromoterFields", "MotifFields", "CoexFields"]


def getDistinctFields(cur, tableName):
	cur.execute('SELECT DISTINCT Field FROM ' + tableName)
	result = []
	for record in cur.fetchall():
		result += record
	return result

def calculateOverlap(cur, table1, field1, table2, field2):
	cur.execute("SELECT COUNT(Gene) FROM " + table1 + " WHERE Field = '" + field1 + "'")
	count1 = int(cur.fetchone()[0])
	cur.execute("SELECT COUNT(Gene) FROM " + table2 + " WHERE Field = '" + field2 + "'")
	count2 = int(cur.fetchone()[0])
	cur.execute("SELECT COUNT(Gene) FROM " + table1 + " WHERE Field = '" + field1 + "' AND Gene IN (SELECT Gene FROM " + table2 + " WHERE Field = '" + field2 + "')")
	commonCount = int(cur.fetchone()[0])
	overlapMin = float(commonCount) / float(max(count1, count2))
	overlapMax = float(commonCount) / float(min(count1, count2))
	return (commonCount, overlapMin, overlapMax)

def overlapTables(cur, tableName1, tableName2):
	fieldNames1 = getDistinctFields(cur, tableName1)
	fieldNames2 = getDistinctFields(cur, tableName2)
	for f1 in fieldNames1:
		for f2 in fieldNames2:
			(commonCount, overlapMin, overlapMax) = calculateOverlap(cur, tableName1, f1, tableName2, f2)
			if overlapMax < 0.5:
				continue
			print(tableName1 + ':' + f1, tableName2 + ':' + f2, commonCount, overlapMin, overlapMax, sep='\t')
			
	
	
print('Field A', 'Field B', 'Common genes', 'Min Overlap', 'Max Overlap', sep='\t')

conn = sqlite3.connect('../Results/yeast.sqlite')

cur = conn.cursor()



for i in range(0, len(tables)):
	j = i+1
	while j < len(tables):
		overlapTables(cur, tables[i], tables[j])
		j = j+1

#conn.commit()
conn.close()