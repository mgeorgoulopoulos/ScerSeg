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
This utility will attempt to patch missing values in given CSV file. See below
for description of the algorithm.
This program expects 2 arguments: inputCsv and outputCsv.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "utils/CsvReader.h"

using namespace std;
using namespace Csv;

// Returns all the rows that have defined the cells that rowIndex has defined.
vector<int> findCompatibleRows(const Table &table, int rowIndex) {
	const Row &row = table[rowIndex];

	vector<int> result;
	for (int i = 1; i < (int)table.size(); i++) {
		if (i == rowIndex) {
			continue;
		}

		const Row &candidateRow = table[i];

		bool accepted = true;
		for (int j = 1; j < (int)row.size(); j++) {
			if (row[j].type == Cell::Number) {
				if (candidateRow[j].type != Cell::Number) {
					accepted = false;
					break;
				}
			}
		}

		if (accepted) {
			result.push_back(i);
		}
	}

	return result;
}

// Patch cell in given row/column using values from compatibleRows.
double patchCell(const Table &table, int rowIndex, int columnIndex,
				 const vector<int> &compatibleRows) {

	const Row &row = table[rowIndex];

	// get the subset of compatible rows which have the missing column defined
	vector<int> compatibleRowsFilled;
	for (const int i : compatibleRows) {
		if (table[i][columnIndex].type != Cell::Number) {
			continue;
		}
		compatibleRowsFilled.push_back(i);
	}

	if (compatibleRowsFilled.empty()) {
		printf("Patching with 0 for %s : %s\n",
			   cellToString(table[rowIndex][0]).c_str(),
			   cellToString(table[0][columnIndex]).c_str());
		return 0.0;
	}

	// Sort the compatible rows using Euclidean distance
	auto lessThanEuclidean = [&](int a, int b) {
		const Row &rowA = table[a];
		const Row &rowB = table[b];

		double distA = 0.0;
		double distB = 0.0;

		for (int i = 1; i < (int)row.size(); i++) {
			if (row[i].type != Cell::Number) {
				continue;
			}
			const double axisDistA = rowA[i].number - row[i].number;
			distA += axisDistA * axisDistA;

			const double axisDistB = rowB[i].number - row[i].number;
			distB += axisDistB * axisDistB;
		}

		return distA < distB;
	};
	std::sort(compatibleRowsFilled.begin(), compatibleRowsFilled.end(),
			  lessThanEuclidean);

	// Keep at most 3 elements
	if ((int)compatibleRowsFilled.size() > 3) {
		compatibleRowsFilled.resize(3);
	}

	// Now sort again using the value of the column
	auto lessThan = [&](int a, int b) {
		return table[a][columnIndex].number < table[b][columnIndex].number;
	};
	std::sort(compatibleRowsFilled.begin(), compatibleRowsFilled.end(),
			  lessThan);
	const int medianRowIndex =
		compatibleRowsFilled[compatibleRowsFilled.size() / 2];

	printf("[%d, %d]: using: ", rowIndex, columnIndex);
	for (int i : compatibleRowsFilled) {
		printf("%d, ", i);
	}
	printf("[%f]\n", table[medianRowIndex][columnIndex].number);

	// Done. Return median
	return table[medianRowIndex][columnIndex].number;
}

// Finds missing values in the table and attempts to patch them
bool patchMissingValues(Table *table) {
	// Keep an unmodified copy
	const Table originalTable = *table;

	for (int rowIndex = 1; rowIndex < (int)originalTable.size(); rowIndex++) {
		const Row &row = originalTable[rowIndex];

		// Get all compatible rows to fill-in the missing data.
		vector<int> compatibleRows =
			findCompatibleRows(originalTable, rowIndex);

		for (int columnIndex = 1; columnIndex < (int)row.size();
			 columnIndex++) {
			if (row[columnIndex].type != Cell::Null) {
				continue;
			}

			// We have an empty cell here - patch it
			table->at(rowIndex)[columnIndex].type = Cell::Number;
			table->at(rowIndex)[columnIndex].number =
				patchCell(originalTable, rowIndex, columnIndex, compatibleRows);
		}
	}

	return true;
}

int main(int argc, char *argv[]) {
	string inputFilename = "PrimarySources/SelectedHistonesPromoter.csv";
	string outputFilename = "Results/SelectedHistonesPromoter-Patched.csv";

	if (argc == 3) {
		inputFilename = argv[1];
		outputFilename = argv[2];
	}

	printf("Reading file: %s\n", inputFilename.c_str());
	Table csvTable;
	if (!readCsv(inputFilename, &csvTable)) {
		printf("Failed to load CSV file: %s\n", inputFilename.c_str());
		return 0;
	}

	printf("%d rows. First row will be ignored as header\n",
		   (int)csvTable.size());

	printf("Patching values... ");
	if (!patchMissingValues(&csvTable)) {
		printf("Failed to patch missing values\n");
		return 0;
	}
	printf("Done\n");

	if (!writeCsv(outputFilename, csvTable)) {
		printf("Failed to write output CSV\n");
		return 0;
	}
	printf("Written patched output to %s\n", outputFilename.c_str());

	printf("Full success\n");
	return 0;
}