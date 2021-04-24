#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "glm/glm.hpp"

using namespace std;
using namespace glm;

struct Record {
	int serialNumber = 0;
	string atomName;
	string chain;
	vec3 position;
};

struct Model {
	vector<Record> records;

	struct Chain {
		string name;
		vector<vec3> vertices;
	};

	vector<Chain> chains;

	bool constructChains();
};

bool readAtomRecord(const string &line, Record *record) {
	if (line.size() < 70) {
		//printf("Line short: %s\n", line.c_str());
		return false;
	}

	if (line.substr(0, 4) != "ATOM") {
		//printf("Line not ATOM: %s\n", line.c_str());
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

bool readRecord(const string &line, Record *record) {
	if (line.empty()) {
		return false;
	}

	if (readAtomRecord(line, record)) {
		return true;
	}

	return false;
}

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

// Encapsulates a map from base index to 3D position
struct BasePositionMap {
	double basesPerSegment = 1.0;
	vector<vec3> vertices;

	vec3 position(int indexBase1) const;
    int totalBases() const {
        return vertices.empty()
                    ? 0
                    : (int)(basesPerSegment * vertices.size() - 1);
    }
};

bool extractBasePositionMaps(const Model &model, vector<BasePositionMap> *basePositionMaps) {
	basePositionMaps->clear();

	// Known Scharomyces chromosome sizes, in bp
    vector<int> chromosomeBaseCounts{
        230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643,
        439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066};

	if (model.chains.size() != chromosomeBaseCounts.size()) {
        printf("Unexpected number of chromosomes: Chains=%d , known "
                "chromosome lengths=%d\n",
                model.chains.size(), chromosomeBaseCounts.size());
    }

	for (int i = 0; i < static_cast<int>(model.chains.size()); i++) {
		const Model::Chain &chain = model.chains[i];
		BasePositionMap m;
		m.vertices = chain.vertices;
		const int bp = chromosomeBaseCounts[i];
		const int segmentCount = static_cast<int>(chain.vertices.size()) - 1;
		if (segmentCount <= 0) {
			printf("W: Chain %s has zero vertices. Ignore this chain.\n", chain.name.c_str());
			basePositionMaps->push_back(m);
			continue;
		}
		m.basesPerSegment = (double)bp / (double)segmentCount;
		basePositionMaps->push_back(m);
	}

	assert(basePositionMaps->size() == model.chains.size());
	
	return true;
}

using Row = vector<string>;
using Table = vector<Row>;

bool readTSV(const string &filename, Table *table) {
	table->clear();

	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == nullptr) {
		printf("Failed to open TSV file: %s\n", filename.c_str());
		return false;
	}

	// get file size
	fseek(fp, 0, SEEK_END);
	const int fileSize = (int) ftell(fp);
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

struct Locus {
	string geneName;
	int chromosomeBase1 = 1;
	int baseStart = 1;
	int baseEnd = 1;
	vec3 position;
};

bool readLoci(const string &filename, vector<Locus> *loci) {
	loci->clear();

	Table table;
	if (!readTSV(filename, &table)) {
		return false;
	}

	loci->reserve(table.size());

	for (int i = 1 /* skip column names */; i < (int)table.size(); i++) {
		const Row &row = table[i];
		if (row.size() != 4) {
			printf("Loci file (%s) line %d has %d values. We only support 4.\n", filename.c_str(), i, row.size());
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
	const string filename = "41586_2010_BFnature08973_MOESM239_ESM.pdb";

	Model model;
	if (!readModel(filename, &model)) {
		printf("Failed to read PDB file: %s\n", filename.c_str());
		return 0;
	}

	printf("Read %d records from PDB file.\n", model.records.size());
	printf("Read %d chains from PDB file.\n", model.chains.size());

	vector<BasePositionMap> basePositionMaps;
	if (!extractBasePositionMaps(model, &basePositionMaps)) {
		printf("Failed to extract base-position maps from 3D model\n");
		return 0;
	}

	// Print some statistics
	vec3 minPos, maxPos, avgPos;
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
		avgPos /= model.records.size();
	}
	printf("Min pos: %.03f\t%.03f\t%.03f\n", minPos.x, minPos.y, minPos.z);
	printf("Max pos: %.03f\t%.03f\t%.03f\n", maxPos.x, maxPos.y, maxPos.z);
	printf("Avg pos: %.03f\t%.03f\t%.03f\n", avgPos.x, avgPos.y, avgPos.z);


    // How far apart is 1 kb in space? Take position samples every 1000 base
    // pairs and output to a file to do some statistics later on.
	{
		const string kbFilename = "kb.csv";
		FILE *fp = fopen(kbFilename.c_str(), "w");
		if (fp == nullptr) {
			printf("Failed to open %s for writing\n", kbFilename.c_str());
			return 0;
		}
		for (const auto &m : basePositionMaps) {
			for (int i = 1; i < m.totalBases() - 1000; i += 1000) {
				vec3 a = m.position(i);
				vec3 b = m.position(i + 1);
				const double distance = glm::distance(a, b);
				fprintf(fp, "%f\n", distance);
			}
		}
		fclose(fp);

		printf("Written kilobase distance measurements to %s\n", kbFilename.c_str());
	}

    // Now, load gene loci (start/end positions in chromosome).
	vector<Locus> loci;
	if (!readLoci("Loci.tsv", &loci)) {
		printf("Failed to load Loci file\n");
		return 0;
	}

	printf("Read loci for %d genes\n", loci.size());
	
	// Finally, map each gene to a position in space. Use median base.
	for (Locus &gene : loci) {
		const int medianBase = (gene.baseStart + gene.baseEnd) / 2;
		const int chromosomeIndex = gene.chromosomeBase1 - 1;
		if (chromosomeIndex < 0 || chromosomeIndex >= (int)basePositionMaps.size()) {
			printf("Invalid chromosome base 1 index: %d (Gene %s)\n", gene.chromosomeBase1, gene.geneName.c_str());
			return 0;
		}
		const BasePositionMap &m = basePositionMaps[chromosomeIndex];
		gene.position = m.position(medianBase);
	}

	// Finally - write out the gene positions
	const string outputFilename = "GenePositions.tsv";
	FILE *fp = fopen(outputFilename.c_str(), "w");
	if (fp == nullptr) {
		printf("Failed to open %s for writing\n", outputFilename.c_str());
		return 0;
	}
	fprintf(fp, "Gene\tChromosome\tStart\tEnd\tx\ty\tz\n");
	for (const Locus &gene : loci) {
		fprintf(fp, "%s\t%d\t%d\t%d\t%.03f\t%.03f\t%.03f\n",
			gene.geneName.c_str(),
			gene.chromosomeBase1, gene.baseStart, gene.baseEnd,
			gene.position.x, gene.position.y, gene.position.z);
	}
	fclose(fp);

	printf("Gene positions written to %s\n", outputFilename.c_str());
	printf("Full success.\n");

    return 0;

}

bool Model::constructChains()
{
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

vec3 BasePositionMap::position(int indexBase1) const {
	if (vertices.empty()) {
		return vec3();
	}

	const int indexSpace = std::max(0, indexBase1 - 1);
	const double segmentSpace = indexSpace / basesPerSegment;
	int vertexA = (int) segmentSpace;
	vertexA = std::min((int)vertices.size() - 1, vertexA);
	int vertexB = vertexA + 1;
	vertexB = std::min((int)vertices.size() - 1, vertexB);
	if (vertexB <= vertexA) {
		return vertices[vertexA];
	}

	// Interpolate to get 3D position of base index
	double tmp;
	const double t = modf(segmentSpace, &tmp);
	return mix(vertices[vertexA], vertices[vertexB], t);

}