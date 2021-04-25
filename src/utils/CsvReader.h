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
This module handles loading/saving of CSV files.
The files are expected to be well-behaved: same number of columns in all rows.
No escape characters etc.
*/

#include <stdio.h>
#include <string>
#include <vector>

namespace Csv {

using namespace std;

bool isNumber(const std::string &s, double *value) {
	double num = 0.0;
	try {
		num = std::stod(s);
	} catch (...) {
		return false;
	}

	*value = num;
	return true;
}

struct Cell {
	enum Type { Null, Text, Number };
	string text;
	double number;
	Type type = Null;
};

Cell stringToCell(const string &cellText) {
	Cell result;
	if (cellText.empty()) {
		result.type = Cell::Null;
		return result;
	}

	double num = 0.0;
	if (isNumber(cellText, &num)) {
		result.number = num;
		result.type = Cell::Number;
		return result;
	}

	result.type = Cell::Text;
	result.text = cellText;
	return result;
}

string cellToString(const Cell &cell) {
	if (cell.type == Cell::Null) {
		return "";
	}

	if (cell.type == Cell::Text) {
		return cell.text;
	}

	char buffer[64];
	sprintf(buffer, "%.9f", cell.number);
	return buffer;
}

using Row = vector<Cell>;
using Table = vector<Row>;

// Minimal error checking here: the CSV is required to be well-behaved: comma
// means comma at all times.
bool readCsv(const string &filename, Table *table) {
	table->clear();

	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == nullptr) {
		printf("Error: Failed to open %s for reading\n", filename.c_str());
		return false;
	}

	fseek(fp, 0, SEEK_END);
	const int fileSize = (int)ftell(fp);
	fseek(fp, 0, SEEK_SET);

	Row row;
	string token;
	bool firstRow = true;
	int rowSize = 0;
	for (int i = 0; i < fileSize; i++) {
		const char c = getc(fp);
		if (c == 0x0A || c == 0x0D) {
			if (!row.empty()) {
				// Technically not correct - this won't catch the empty line
				// denoting a single empty row, but, oh well...
				row.push_back(stringToCell(token));
				token.clear();
			}
			if (row.empty()) {
				continue;
			}
			if (firstRow) {
				firstRow = false;
				rowSize = (int)row.size();
			} else {
				if ((int)row.size() != rowSize) {
					printf("Error: Found row of size %d which is not equal to "
						   "first row "
						   "(%d)\n",
						   (int)row.size(), rowSize);
					fclose(fp);
					return false;
				}
			}
			table->push_back(row);
			row.resize(0);
			continue;
		} // end if (newline)

		if (c == ',') {
			row.push_back(stringToCell(token));
			token.clear();
			continue;
		}

		token += c;
	}
	if (!token.empty()) {
		row.push_back(stringToCell(token));
	}
	if (!row.empty()) {
		if (firstRow) {
			firstRow = false;
			rowSize = (int)row.size();
		} else {
			if ((int)row.size() != rowSize) {
				printf("Error: Found row of size %d which is not equal to "
					   "first row "
					   "(%d)\n",
					   (int)row.size(), rowSize);
				fclose(fp);
				return false;
			}
		}
		table->push_back(row);
		row.resize(0);
	}

	fclose(fp);

	return true;
}

bool writeCsv(const string &filename, const Table &table) {
	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == nullptr) {
		printf("Error: Failed to open %s for writing\n", filename.c_str());
		return false;
	}
	for (const Row &row : table) {
		if (row.empty()) {
			fprintf(fp, "\n");
			continue;
		}
		fprintf(fp, "%s", cellToString(row[0]).c_str());
		for (int i = 1; i < (int)row.size(); i++) {
			const Cell &cell = row[i];
			fprintf(fp, ",%s", cellToString(cell).c_str());
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
}

} // end namespace Csv