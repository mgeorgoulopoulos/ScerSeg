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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "utils/TsvReader.h"
#include "utils/Vec3D.h"

using namespace std;

// Represents a record of a PDB file. We only keep the fields relevant to yeast
// 3D model.
struct Record {
	int serialNumber = 0;
	string atomName;
	string chain;
	Vec3D position;
};

// Represents the yeast 3D model. We break the records in chains (which are
// chromosomes) and we also keep all the records in a central container.
struct Model {
	// All the records (definition points of the 3D model).
	vector<Record> records;

	struct Chain {
		// Chromosome name
		string name;

		// Definition points
		vector<Vec3D> vertices;
	};

	// All the chromosomes, each represented by a list of definition points.
	vector<Chain> chains;

	bool constructChains();
};

// Parses one line of PDB file text to a record
bool readAtomRecord(const string &line, Record *record) {
	if (line.size() < 70) {
		// printf("Line short: %s\n", line.c_str());
		return false;
	}

	if (line.substr(0, 4) != "ATOM") {
		// printf("Line not ATOM: %s\n", line.c_str());
		return false;
	}

	const auto serialNumberString = line.substr(6, 5);
	record->serialNumber = atoi(serialNumberString.c_str());

	const auto atomName = line.substr(12, 4);
	record->atomName = "";
	for (const char c : atomName) {
		if (c != ' ') {
			record->atomName += c;
		}
	}

	record->chain = "";
	record->chain += line[21];

	const auto xString = line.substr(30, 8);
	record->position.x = atof(xString.c_str());

	const auto yString = line.substr(38, 8);
	record->position.y = atof(yString.c_str());

	const auto zString = line.substr(46, 8);
	record->position.z = atof(zString.c_str());

	return true;
}

// Reads one line from the file and parses the line to a record
bool readRecord(const string &line, Record *record) {
	if (line.empty()) {
		return false;
	}

	if (readAtomRecord(line, record)) {
		return true;
	}

	return false;
}

// Read the complete 3D model from the PDB file.
bool readModel(const string &filename, Model *model) {
	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == nullptr) {
		printf("Failed to open PDB file: %s\n", filename.c_str());
		return false;
	}

	// Get file size
	fseek(fp, 0, SEEK_END);
	const int fileSize = static_cast<int>(ftell(fp));
	fseek(fp, 0, SEEK_SET);

	// Read 70-character records
	string line;
	for (int i = 0; i < fileSize; i++) {
		const char c = getc(fp);
		if (c == 0x0A || c == 0x0D) {
			Record r;
			if (readRecord(line, &r)) {
				model->records.push_back(r);
			}
			line.clear();
			continue;
		}
		line.push_back(c);
	}

	{ // handle last line
		Record r;
		if (readRecord(line, &r)) {
			model->records.push_back(r);
		}
		line.clear();
	}

	fclose(fp);

	if (!model->constructChains()) {
		return false;
	}

	return true;
}

// Encapsulates a map from base-pair index to 3D position. Base-pair indices
// index one chromosome and their domain is [1,chromosomeSizeInBasePairs].
struct BasePositionMap {
	double basesPerSegment = 1.0;
	vector<Vec3D> vertices;

	// Calculates interpolated position of base-pair index.
	Vec3D position(int indexBase1) const;

	int totalBases() const {
		return vertices.empty() ? 0
								: (int)(basesPerSegment * vertices.size() - 1);
	}
};

// Converts a Model to a list of BasePositionMap objects. For yeast we will get
// a list of 16 BasePositionMap objects. So when we want to get the 3D position
// of chromosome i, basepair j: we will query the i-th BasePositionMap for its
// j-th base-pair position. ASSUMPTION: We assume that the control points in the
// PDB file are evenly distributed in the base-pair domain.
bool extractBasePositionMaps(const Model &model,
							 vector<BasePositionMap> *basePositionMaps) {
	basePositionMaps->clear();

	// Known Scharomyces chromosome sizes, in bp
	vector<int> chromosomeBaseCounts{
		230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643,
		439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066};

	if (model.chains.size() != chromosomeBaseCounts.size()) {
		printf("Unexpected number of chromosomes: Chains=%d , known "
			   "chromosome lengths=%d\n",
			   (int)model.chains.size(), (int)chromosomeBaseCounts.size());
	}

	for (int i = 0; i < static_cast<int>(model.chains.size()); i++) {
		const Model::Chain &chain = model.chains[i];
		BasePositionMap m;
		m.vertices = chain.vertices;
		const int bp = chromosomeBaseCounts[i];
		const int segmentCount = static_cast<int>(chain.vertices.size()) - 1;
		if (segmentCount <= 0) {
			printf("W: Chain %s has zero vertices. Ignore this chain.\n",
				   chain.name.c_str());
			basePositionMaps->push_back(m);
			continue;
		}
		m.basesPerSegment = (double)bp / (double)segmentCount;
		basePositionMaps->push_back(m);
	}

	assert(basePositionMaps->size() == model.chains.size());

	return true;
}

// A continuous area in the genome. Located at a specific chromosome and
// specifying its start and end base pair index. All indices are 1-based.
struct Locus {
	string geneName;

	// Index of chromosome. 1 = ChrI
	int chromosomeBase1 = 1;

	// Start and end base-pair index
	int baseStart = 1;
	int baseEnd = 1;

	// Position in space. This refers to the interpolated position based on the
	// yeast 3D model and computed for the center base-pair of the gene. So if a
	// gene is 1000 bases long, we calculate the 3D position of base 500.
	Vec3D position;
};

// Reads tab-separated file with columns: gene name, chromosome, start base, end
// base. Returns a list of Locus objects with an uninitialized 3D position. This
// will be filled later on based on their chromosome and base-pair indices.
bool readLoci(const string &filename, vector<Locus> *loci) {
	using namespace Tsv;

	loci->clear();

	Table table;
	if (!readTSV(filename, &table)) {
		return false;
	}

	loci->reserve(table.size());

	for (int i = 1 /* skip column names */; i < (int)table.size(); i++) {
		const Row &row = table[i];
		if (row.size() != 4) {
			printf("Loci file (%s) line %d has %d values. We only support 4.\n",
				   filename.c_str(), i, (int)row.size());
			return false;
		}

		Locus locus;
		locus.geneName = row[0];
		locus.chromosomeBase1 = atoi(row[1].c_str());
		if (locus.chromosomeBase1 <= 0) {
			printf("Invalid chromosome number %s : %d\n", filename.c_str(), i);
			return false;
		}
		locus.baseStart = atoi(row[2].c_str());
		if (locus.baseStart <= 0) {
			printf("Invalid start base %s : %d\n", filename.c_str(), i);
			return false;
		}
		locus.baseEnd = atoi(row[3].c_str());
		if (locus.baseEnd <= 0) {
			printf("Invalid base end %s : %d\n", filename.c_str(), i);
			return false;
		}

		loci->push_back(locus);
	}

	return true;
}

int main(int argc, char *argv[]) {
	// Load the 3D model - I'm not adding it in the package as it is not my own
	// work. Notifying the user that she must find and add it.
	Model model;
	const string filename =
		"PrimarySources/41586_2010_BFnature08973_MOESM239_ESM.pdb";
	if (!readModel(filename, &model)) {
		printf("Failed to read PDB file: %s. If you haven't done so already, "
			   "please find this file from the corresponding paper by Duan et "
			   "al. and place it in 'PrimarySources' folder of this package.\n",
			   filename.c_str());
		return 0;
	}

	printf("Read %d records from PDB file.\n", (int)model.records.size());
	printf("Read %d chains from PDB file.\n", (int)model.chains.size());

	// Convert the model to list of BasePositionMap obejcts so we can map our
	// gene positions.
	vector<BasePositionMap> basePositionMaps;
	if (!extractBasePositionMaps(model, &basePositionMaps)) {
		printf("Failed to extract base-position maps from 3D model\n");
		return 0;
	}

	// Print some statistics
	Vec3D minPos, maxPos, avgPos;
	bool first = true;
	for (const auto &r : model.records) {
		if (first) {
			first = false;
			minPos = maxPos = avgPos = r.position;
			continue;
		}
		minPos.x = std::min(minPos.x, r.position.x);
		minPos.y = std::min(minPos.y, r.position.y);
		minPos.z = std::min(minPos.z, r.position.z);
		maxPos.x = std::max(maxPos.x, r.position.x);
		maxPos.y = std::max(maxPos.y, r.position.y);
		maxPos.z = std::max(maxPos.z, r.position.z);
		avgPos += r.position;
	}
	if (model.records.size() > 0) {
		avgPos /= (double)model.records.size();
	}
	printf("Min pos: %.03f\t%.03f\t%.03f\n", minPos.x, minPos.y, minPos.z);
	printf("Max pos: %.03f\t%.03f\t%.03f\n", maxPos.x, maxPos.y, maxPos.z);
	printf("Avg pos: %.03f\t%.03f\t%.03f\n", avgPos.x, avgPos.y, avgPos.z);

#if 0 // Set this to 1 to write kb.csv containing sampled kilobase distances in
	  // 3D space. This might be useful to get an idea of the units the model is
	  // in. They are not in any standard units as far as I know.
	// How far apart is 1 kb in space? Take position samples every 1000 base
	// pairs and output to a file to do some statistics later on.
	{
		const string kbFilename = "Results/kb.csv";
		FILE *fp = fopen(kbFilename.c_str(), "w");
		if (fp == nullptr) {
			printf("Failed to open %s for writing\n", kbFilename.c_str());
			return 0;
		}
		for (const auto &m : basePositionMaps) {
			for (int i = 1; i < m.totalBases() - 1000; i += 1000) {
				Vec3D a = m.position(i);
				Vec3D b = m.position(i + 1);
				const double distance = Vec3D::distance(a, b);
				fprintf(fp, "%f\n", distance);
			}
		}
		fclose(fp);

		printf("Written kilobase distance measurements to %s\n",
			   kbFilename.c_str());
	}
#endif // Write kilobase distances file

	// Now, load gene loci (start/end positions in chromosome).
	vector<Locus> loci;
	if (!readLoci("PrimarySources/Loci.tsv", &loci)) {
		printf("Failed to load Loci file\n");
		return 0;
	}

	printf("Read loci for %d genes\n", (int)loci.size());

	// Finally, map each gene to a position in space. Use median base.
	for (Locus &gene : loci) {
		const int medianBase = (gene.baseStart + gene.baseEnd) / 2;
		const int chromosomeIndex = gene.chromosomeBase1 - 1;
		if (chromosomeIndex < 0 ||
			chromosomeIndex >= (int)basePositionMaps.size()) {
			printf("Invalid chromosome base 1 index: %d (Gene %s)\n",
				   gene.chromosomeBase1, gene.geneName.c_str());
			return 0;
		}
		const BasePositionMap &m = basePositionMaps[chromosomeIndex];
		gene.position = m.position(medianBase);
	}

	// Finally - write out the gene positions
	const string outputFilename = "Results/GenePositions.tsv";
	FILE *fp = fopen(outputFilename.c_str(), "w");
	if (fp == nullptr) {
		printf("Failed to open %s for writing\n", outputFilename.c_str());
		return 0;
	}
	fprintf(fp, "Gene\tChromosome\tStart\tEnd\tx\ty\tz\n");
	for (const Locus &gene : loci) {
		fprintf(fp, "%s\t%d\t%d\t%d\t%.03f\t%.03f\t%.03f\n",
				gene.geneName.c_str(), gene.chromosomeBase1, gene.baseStart,
				gene.baseEnd, gene.position.x, gene.position.y,
				gene.position.z);
	}
	fclose(fp);

	printf("Gene positions written to %s\n", outputFilename.c_str());
	printf("Full success.\n");

	return 0;
}

bool Model::constructChains() {
	chains.clear();

	map<string, Chain> nameToChain;
	for (const Record &r : records) {
		nameToChain[r.chain].vertices.push_back(r.position);
	}

	for (const auto it : nameToChain) {
		Chain chain = it.second;
		chain.name = it.first;
		chains.push_back(chain);
	}

	return true;
}

// This is the important part: mapping base indices to positions. We do so by
// mapping the base index between two control points of the original model,
// producing a parameter t that tells us exactly how distant the base pair is
// from the first and the second control point (0 = on top of the first point, 1
// = on the second point, 0.3 = 30% of the distance between control points).
Vec3D BasePositionMap::position(int indexBase1) const {
	if (vertices.empty()) {
		return Vec3D();
	}

	const int indexSpace = std::max(0, indexBase1 - 1);
	const double segmentSpace = indexSpace / basesPerSegment;
	int vertexA = (int)segmentSpace;
	vertexA = std::min((int)vertices.size() - 1, vertexA);
	int vertexB = vertexA + 1;
	vertexB = std::min((int)vertices.size() - 1, vertexB);
	if (vertexB <= vertexA) {
		return vertices[vertexA];
	}

	// Interpolate to get 3D position of base index
	double tmp;
	const double t = modf(segmentSpace, &tmp);
	return Vec3D::mix(vertices[vertexA], vertices[vertexB], t);
}