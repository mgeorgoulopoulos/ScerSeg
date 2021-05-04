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
This module handles loading of TSV files.
The files are expected to be well-behaved: same number of columns in all rows.
No escape characters etc.
*/

#include <stdio.h>
#include <string>
#include <vector>

namespace Tsv {

using namespace std;
using Row = vector<string>;
using Table = vector<Row>;

bool writeTsv(const string &filename, const Table &table) {
	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == nullptr)
		return false;
	for (const Row &row : table) {
		for (int i = 0; i < (int)row.size() - 1; i++) {
			fprintf(fp, "%s\t", row[i].c_str());
		}
		if (!row.empty())
			fprintf(fp, row.back().c_str());
		fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
}

bool readTSV(const string &filename, Table *table) {
	table->clear();

	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == nullptr) {
		printf("Failed to open TSV file: %s\n", filename.c_str());
		return false;
	}

	// get file size
	fseek(fp, 0, SEEK_END);
	const int fileSize = (int)ftell(fp);
	fseek(fp, 0, SEEK_SET);

	string token;
	Row row;
	for (int i = 0; i < fileSize; i++) {
		const char c = getc(fp);
		if (c == 0x0A || c == 0x0D) {
			if (!token.empty()) {
				row.push_back(token);
			}
			if (!row.empty()) {
				table->push_back(row);
				token.clear();
				row.resize(0);
			}
			continue;
		}
		if (c == '\t') {
			row.push_back(token);
			token.clear();
			continue;
		}
		token += c;
	}
	if (!row.empty()) {
		table->push_back(row);
		row.clear();
	}

	fclose(fp);

	return true;
}

} // end namespace Tsv