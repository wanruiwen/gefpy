
#include "CellExpWriter.h"

#include <iostream>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>



CellExpWriter::CellExpWriter(const string &inPath, const string &outPath) {
    bool r = copyFile(inPath, outPath);
    if(!r) cerr << "Error writing to file (" << outPath << ") is failed!" << endl;

    printf("create h5 file: %s\n", inPath.c_str());
    m_file_id = H5Fcreate(outPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT);

//    storeVersion();
    m_group_id = H5Gcreate(m_file_id, "/cellExp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_status = H5Gclose(m_group_id);
}

CellExpWriter::~CellExpWriter()
{
    m_status = H5Fclose(m_file_id);
}

bool CellExpWriter::storeVersion() {
    hsize_t dimsAttr[1] = {1};
    hid_t data_space = H5Screate_simple(1, dimsAttr, NULL);
    hid_t attr = H5Acreate(m_file_id, "version", H5T_STD_U32LE, data_space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_UINT, &version);
    H5Sclose(data_space);
    H5Aclose(attr);
    return true;
}

bool CellExpWriter::storeGeneList(vector<string> & geneList) {
    hsize_t dims[1] = {geneList.size()};

    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 32);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(m_file_id, "cellExp/geneList", strtype, dataspace_id, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    m_status = H5Dwrite(dataset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &geneList[0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return true;
}


bool CellExpWriter::copyFile(const string &inPath, const string &outPath) {
    ifstream fin(inPath, ios::binary);
    ofstream fout(outPath, ios::binary);

    bool bRet = true;

    while(!fin.eof()){
        char szBuf;
        fin.read(&szBuf, sizeof(char));

        if(fin.eof()) break;

        if (fout.bad())
        {
            bRet = false;
            break;
        }
        fout.write(&szBuf, sizeof(char));
    }

    fout.close();
    fin.close();
    return bRet;
}

bool CellExpWriter::storeCellBorder(byte borderPath[][16][2], int size) {
    hsize_t dims[3];
    dims[0] = size;
    dims[1] = 16;
    dims[2] = 2;

    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
    hid_t dataset_id = H5Dcreate(m_file_id, "cellExp/cellBorder", H5T_STD_U8LE, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &borderPath[0][0][0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return true;
}

static bool myCompare(const CellExpData& a1,const CellExpData& a2)
{
    return a1.geneID <= a2.geneID;
}


bool CellExpWriter::storeCellExp() {

    hsize_t dims[1] = {cell_exp_list.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, geneID), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);

    filetype = H5Tcreate(H5T_COMPOUND, 4);
    H5Tinsert(filetype, "geneID", 0, H5T_STD_U16LE);
    H5Tinsert(filetype, "count", 2, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(m_file_id, "cellExp/cellExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_exp_list[0]);

    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return true;
}

void CellExpWriter::storeCell(unsigned int * x, unsigned int * y, unsigned short * area, int size) {
    vector<CellData> cellDataList(size);

    for(int i=0; i < size; i++){
        CellData cell = {x[i], y[i], area[i], cell_gene_count_list[i],
                         cell_exp_count_list[i], cell_gene_exp_list[i]};
        cellDataList.emplace_back(cell);
    }

    hsize_t dims[1] = {(hsize_t)size};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellData));
    H5Tinsert(memtype, "x", HOFFSET(CellData, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(CellData, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "area", HOFFSET(CellData, area), H5T_NATIVE_SHORT);
    H5Tinsert(memtype, "geneCount", HOFFSET(CellData, geneCount), H5T_NATIVE_SHORT);
    H5Tinsert(memtype, "expCount", HOFFSET(CellData, expCount), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "geneExp", HOFFSET(CellData, geneExp), H5T_NATIVE_UINT);

    filetype = H5Tcreate(H5T_COMPOUND, 18);
    H5Tinsert(filetype, "x", 0, H5T_STD_U32LE);
    H5Tinsert(filetype, "y", 4, H5T_STD_U32LE);
    H5Tinsert(filetype, "area", 8, H5T_STD_U16LE);
    H5Tinsert(filetype, "geneCount", 10, H5T_STD_U16LE);
    H5Tinsert(filetype, "expCount", 12, H5T_STD_U16LE);
    H5Tinsert(filetype, "geneExp", 14, H5T_STD_U32LE);


    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(m_file_id, "cellExp/cellExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cellDataList[0]);

    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
//    hdset_reg_ref_t wdata[1];
//    status = H5Sselect_elements (space, H5S_SELECT_SET, 4, coords[0]);
//    status = H5Rcreate (&wdata[0], file, DATASET2, H5R_DATASET_REGION, space);
}


void CellExpWriter::cell_bin(unsigned int ** cell_exp_index, unsigned int size){
    unsigned long long xy_bin;
    vector<CellExpData> mergeCellExpTmp;

    for (unsigned int i = 0; i < size; ++i) {
        xy_bin = cell_exp_index[i][0];
        xy_bin = xy_bin << 32 | cell_exp_index[i][1];
        auto iter = gene_exp_map.find(xy_bin);
        if(iter != gene_exp_map.end()){
            vector<CellExpData> cxp = iter->second;
            mergeCellExpTmp.insert(mergeCellExpTmp.end(),cxp.begin(),cxp.end());
        }
    }

    sort(mergeCellExpTmp.begin(), mergeCellExpTmp.end(),myCompare);

    unsigned short expCount = 0;
    vector<CellExpData> cellExp;
    auto it = mergeCellExpTmp.begin();
    while (it != mergeCellExpTmp.end()) {
        if(!cellExp.empty() && it->geneID == cellExp.back().geneID){
            cellExp.back().count += it->count;
            expCount += it->count;
        }else{
            cellExp.emplace_back(*it);
            expCount += it->count;
        }
        ++it;
    }

    unsigned short geneCount = cellExp.size();

    cell_exp_list.insert(cell_exp_list.end(), cellExp.begin(), cellExp.end());
    cell_exp_count_list.emplace_back(expCount);
    cell_gene_count_list.emplace_back(geneCount);
    if(cell_gene_exp_list.empty()){
        cell_gene_exp_list.emplace_back(0);
    }else {
        cell_gene_exp_list.emplace_back(cell_gene_exp_list.back()+geneCount);
    }

    vector<CellExpData>().swap(mergeCellExpTmp);
    vector<CellExpData>().swap(cellExp);
}

void CellExpWriter::setGeneExpMap(const string &inPath){
    GeneExp gene_exp = GeneExp(inPath, 1);
    auto * gene_index = (unsigned int *)malloc(gene_exp.getExpLen() * sizeof(unsigned int));
    vector<string> gene_names;
    gene_exp.getGeneData(gene_index, gene_names);

    Gene * expData = gene_exp.getExpData();

    unsigned long long xy_bin;
    for (int i = 0; i < gene_exp.getExpLen(); ++i) {
        CellExpData cellExpData = {(unsigned short)gene_index[i], (unsigned short)expData[i].cnt};
        xy_bin = expData[i].x;
        xy_bin = xy_bin << 32 | expData[i].y;

        auto iter = gene_exp_map.find(xy_bin);
        if(iter != gene_exp_map.end()){
            iter->second.emplace_back(cellExpData);
        }else{
            vector<CellExpData> cellExpDataList;
            cellExpDataList.emplace_back(cellExpData);
            gene_exp_map.insert(map<unsigned long long, vector<CellExpData>>::value_type (xy_bin, cellExpDataList));
        }
    }

    free(expData);
    free(gene_index);
}





