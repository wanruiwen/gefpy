
#pragma once

#include "hdf5.h"
#include "GeneExp.h"

typedef uint8_t byte;

typedef struct
{
    unsigned int x;
    unsigned int y;
    unsigned short area;
    unsigned short geneCount;
    unsigned int expCount;
    unsigned int geneExp;
//    hdset_reg_ref_t geneExp;
} CellData;

typedef struct
{
    unsigned short geneID;
    unsigned short count;
} CellExpData;

class CellExpWriter
{
public:
    CellExpWriter(const string& inPath, const string& outPath);
    ~CellExpWriter();

    void storeCell(unsigned int * x, unsigned int * y, unsigned short * area, int size);

    bool storeCellExp();

    bool storeCellBorder(byte borderPath[][16][2], int size);

    bool storeGeneList(vector<string>& geneList);

    bool storeVersion();

    static bool copyFile(const string& inPath, const string& outPath);

    void cell_bin(unsigned int ** cell_exp_index, unsigned int size);

    void setGeneExpMap(const string &inPath);

private:
    hid_t m_file_id;
    herr_t m_status;
    hid_t m_group_id;
    map<unsigned long long , vector<CellExpData>> gene_exp_map;
    vector<CellExpData> cell_exp_list;
    vector<unsigned short> cell_gene_count_list;
    vector<unsigned int> cell_exp_count_list;  //offset
    vector<unsigned int> cell_gene_exp_list;  //offset
};