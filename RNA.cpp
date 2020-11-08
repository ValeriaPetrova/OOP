#include <iostream>
#include <stdint.h>
#include <cstring>
#include <cmath>
#include "RNA.h"


RNA::RNA(Nucl N, size_t length) {//+
    for (size_t i = 0; i < length; i++){
        add_elem(N);
    }
}

RNA::RNA(const RNA &rna) {//+
    numb_of_nucl = rna.numb_of_nucl;
    length_of_chain = rna.length_of_chain;
    chain_of_nucl = new size_t [length_of_chain];
    memcpy(chain_of_nucl, rna.chain_of_nucl, sizeof(size_t) * length_of_chain);
}

RNA::~RNA() {//+
    if (chain_of_nucl != nullptr) {
        delete[] chain_of_nucl;
        chain_of_nucl = nullptr;
    }
}

Nucl RNA::GetNucl(size_t idx) const{ //+
    size_t bitMask = 3;
    size_t idxArray = idx / nucl_count_in_size_t;
    size_t idxRNA = idx % nucl_count_in_size_t + 1;
    size_t nucl = (chain_of_nucl[idxArray] & bitMask << (sizeof(size_t) * 8 - idxRNA * 2)) >> (sizeof(size_t) * 8 - idxRNA * 2);
    switch (nucl){
        case 0: return A;
            break;
        case 1: return G;
            break;
        case 2: return C;
            break;
        case 3: return T;
            break;
        default:
            break;
    }
}

void RNA::add_elem(Nucl nucleotide) {//++
    if(chain_of_nucl == nullptr){
        chain_of_nucl = new size_t;
        int shift = 2 * (nucl_count_in_size_t * (length_of_chain + 1) - numb_of_nucl) - 2;
        chain_of_nucl[0] = (size_t)nucleotide << shift;
        length_of_chain++;
    }
    else if (nucl_count_in_size_t > numb_of_nucl / length_of_chain){
        int shift = 2 * (nucl_count_in_size_t * length_of_chain - numb_of_nucl) - 2;
        size_t n = (size_t)nucleotide << shift;
        chain_of_nucl[length_of_chain - 1] = chain_of_nucl[length_of_chain - 1] | n;
        length_of_chain++;
    }
    else { // nucl_count_in_size_t <= numb_of_nucl / length_of_chain
        size_t *new_chain = new size_t[length_of_chain * 2];
        for (int i = 0; i < length_of_chain; i++){
            new_chain[i] = chain_of_nucl[i];
        }
        int shift = 2 * (nucl_count_in_size_t * (length_of_chain + 1) - numb_of_nucl) - 2;
        new_chain[length_of_chain] = (size_t)nucleotide << shift;
        delete[] chain_of_nucl;
        chain_of_nucl = new_chain;
        length_of_chain *= 2;
    }
    numb_of_nucl++;
}

Nucl RNA::complementary(Nucl nucleotide) const{//++
    Nucl n = (Nucl)(3 - nucleotide);
    return n;
}

bool RNA::isComplementary(const RNA& rna) const {//+
    if (this->numb_of_nucl != rna.numb_of_nucl){
        return false;
    } else {
        size_t i = 0;
        while(i != this->numb_of_nucl){
            if (this->GetNucl(i) == !complementary(rna.GetNucl(i))){
                i++;
            } else {
                return false;
            }
        }
        return true;
    }
}

void RNA::trim(size_t last_idx) {  //+ забыть содержимое от lastIndex и дальше
    size_t *new_chain;
    numb_of_nucl = last_idx;
    length_of_chain = last_idx * 2 / sizeof(size_t) / 8 + 1;
    new_chain = chain_of_nucl;
    chain_of_nucl = new size_t[length_of_chain];
    memcpy(chain_of_nucl, new_chain, length_of_chain * sizeof(size_t)); // (new_chain в chain_of_nucl)
    delete[] new_chain;
}

RNA RNA::split(size_t idx) {//+
    RNA result;
    for (size_t i = idx; i < numb_of_nucl; i++) {
        Nucl n = (*this)[i];
        result.add_elem(n);
    }
    trim(idx);
    return result;
}

size_t RNA::cardinality(Nucl value) {
    size_t result = 0;
    for (size_t i = 0; i < this->numb_of_nucl; i++) {
        if ((*this)[i] == value) {
            result++;
        }
    }
    return result;
}

//---------------------------------------------------------------------------------
RNA::reference::reference(size_t idx, RNA &rna1) : num(idx), rna(rna1){ }

RNA operator+(RNA &rna1, RNA &rna2) {//+
    RNA rna;
    rna.numb_of_nucl = rna1.numb_of_nucl + rna2.numb_of_nucl;
    rna.length_of_chain = rna.numb_of_nucl * 2 / sizeof(size_t) / 8 + 1;
    rna.chain_of_nucl = new size_t[rna.length_of_chain];
    memcpy(rna.chain_of_nucl, rna1.chain_of_nucl, rna1.length_of_chain * sizeof(size_t));
    for (size_t i = rna1.numb_of_nucl; i < rna.numb_of_nucl; i++) {
        rna.add_elem(rna2[i]);
    }
    return rna;
}

bool RNA::operator==(const RNA &rna) const{//+
    if (this->numb_of_nucl != rna.numb_of_nucl) {
        return false;
    }
    for (size_t i = 0; i < rna.length_of_chain; i++) {
        if (this->chain_of_nucl[i] != rna.chain_of_nucl[i]) {
            return false;
        }
    }
    return true;
}

bool RNA::operator!=(const RNA &rna) const{//+
    return !((*this) == rna);
}

RNA RNA::operator!(){//+
    for (size_t i = 0; i < this->numb_of_nucl; i++) {
        Nucl nucl = (*this)[i];
        Nucl new_nucl = complementary(nucl);
        (*this)[i] = new_nucl;
    }
    return (*this);
}

RNA& RNA::operator=(RNA const &rna) {//+
    if (this->length_of_chain != rna.length_of_chain) {
        delete[] this->chain_of_nucl;
        this->chain_of_nucl = new size_t[rna.length_of_chain];
        this->length_of_chain = rna.length_of_chain;
    }

    this->numb_of_nucl = rna.numb_of_nucl;
    for (size_t i = 0; i < this->length_of_chain; i++) {
        this->chain_of_nucl[i] = rna.chain_of_nucl[i];
    }
    return (*this);
}

RNA::reference::operator Nucl() const{//+
    return rna.GetNucl(num);
}

RNA::reference RNA::operator[](size_t num) {//+
    return reference(num, *this);
}

RNA::reference& RNA::reference::operator=(Nucl N){//+?
    size_t bitMask = 3;
    size_t idxRNA =  rna.length_of_chain % rna.nucl_count_in_size_t + 1;
    rna.chain_of_nucl[rna.numb_of_nucl] &= ~(bitMask << (sizeof(size_t) * 8 - idxRNA * 2));
    bitMask = N;
    rna.chain_of_nucl[rna.numb_of_nucl] |= (bitMask << (sizeof(size_t) * 8 - idxRNA * 2));
    return *this;
}

RNA::reference& RNA::reference::operator=(const RNA::reference& ref) {//+
    operator=(Nucl(ref));
    return *this;
}

RNA::reference::~reference() = default;

//_________________________________________________________________________________

void RNA::output() const {
    int bits = sizeof(size_t) * 8 - 2;
    size_t bitMask = 3 << bits;
    size_t end_of_rna = 0;
    for (size_t i = 0; i < length_of_chain; i++){
        if (numb_of_nucl - i * nucl_count_in_size_t < nucl_count_in_size_t){
            end_of_rna = numb_of_nucl - i * nucl_count_in_size_t;
        }
        else {
            end_of_rna = nucl_count_in_size_t;
        }
        for (size_t j = 0; j < end_of_rna; j++){
            int n = (int)((chain_of_nucl[i] & (bitMask >> (2 * j))) >> (bits - 2 * j));
            switch (n){
                case 0:
                    std::cout << "A";
                    break;
                case 1:
                    std::cout << "G";
                    break;
                case 2:
                    std::cout << "C";
                    break;
                case 3:
                    std::cout << "T";
                    break;
                default:
                    break;
            }
        }
        std::cout << " ";
    }
}

size_t RNA::length() const {
    return this->length_of_chain;
}

size_t RNA::num_of_nucls() const {
    return this->numb_of_nucl;
}