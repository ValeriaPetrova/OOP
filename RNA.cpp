#include "RNA.h"

RNA::RNA() {}

RNA::RNA(Nucl nucl, size_t length){
    for (size_t i = 0; i < length; i++){
        add_elem(nucl);
    }
}

RNA::RNA(const RNA &rna1) {
    length_of_chain = rna1.length_of_chain;
    chain_of_nucl = new size_t [length_of_chain];
    for (size_t i = 0; i < length_of_chain; i++){
        chain_of_nucl[i] = rna1.chain_of_nucl[i];
    }
    numb_of_nucl = rna1.numb_of_nucl;
}

RNA::~RNA() {
    if (chain_of_nucl != nullptr) {
        delete[] chain_of_nucl;
        chain_of_nucl = nullptr;
    }
}

RNA::reference::reference(size_t idx, RNA &rna1) : num(idx), rna(rna1) {}

RNA::reference::~reference() = default;

//------------------------------------------------------------------------

Nucl RNA::GetNucl(size_t idx) const{
    size_t bitMask = 3;
    size_t idxArray = (size_t)ceil((float)idx / (float)nucl_count_in_size_t) - 1;
    size_t shift = (nucl_count_in_size_t - (idx - nucl_count_in_size_t * idxArray)) * 2;
    size_t nucl = ((*this).chain_of_nucl[idxArray] & (bitMask << shift)) >> shift;
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
    return A;
}

void RNA::add_elem(Nucl nucl){
    if (chain_of_nucl == nullptr){
        chain_of_nucl = new size_t[1];
        length_of_chain++;
        size_t shift = nucl_count_in_size_t * length_of_chain - numb_of_nucl - 1;
        chain_of_nucl[0] = (size_t)nucl << (shift * 2);
    }
    else if (nucl_count_in_size_t > (numb_of_nucl / length_of_chain)){
        size_t shift = nucl_count_in_size_t * length_of_chain - numb_of_nucl - 1;
        size_t n = (size_t)nucl << (shift * 2);
        chain_of_nucl[length_of_chain - 1] = chain_of_nucl[length_of_chain - 1] | n;
    }
    else {
        auto *new_chain = new size_t[length_of_chain + 1];
        for (int i = 0; i < length_of_chain; i++){
            new_chain[i] = chain_of_nucl[i];
        }
        size_t shift = nucl_count_in_size_t * (length_of_chain + 1) - numb_of_nucl - 1;
        new_chain[length_of_chain] = (size_t)nucl << (shift * 2);
        delete[] chain_of_nucl;
        chain_of_nucl = new_chain;
        length_of_chain++;
    }
    numb_of_nucl++;
}

Nucl RNA::complementary(Nucl nucleotide) const{//++
    Nucl n = (Nucl)(3 - nucleotide);
    return n;
}

bool RNA::isComplementary(RNA &rna1) const{
    if (this->num_of_nucls() != rna1.num_of_nucls()) {
        return false;
    }
    RNA buff(rna1);
    return (!buff == (*this));//(!buff == (*this)) ? true : false
}

RNA& RNA::trim(size_t last_idx) {  //+ забыть содержимое от lastIndex и дальше
    size_t *new_chain;
    this->numb_of_nucl = last_idx;
    this->length_of_chain = last_idx * 2 / sizeof(size_t) / 8 + 1;
    new_chain = this->chain_of_nucl;
    this->chain_of_nucl = new size_t[this->length_of_chain];
    memcpy(this->chain_of_nucl, new_chain, this->length_of_chain * sizeof(size_t)); // (new_chain в chain_of_nucl)
    delete[] new_chain;
    return (*this);
}

RNA& RNA::split(size_t idx) {
    RNA rna1;
    rna1.numb_of_nucl = this->numb_of_nucl - idx;
    rna1.length_of_chain = rna1.numb_of_nucl * 2 / sizeof(size_t) / 8 + 1;
    rna1.chain_of_nucl = new size_t[rna1.length_of_chain];
    for (size_t i = idx; i < this->numb_of_nucl; i++) {
        (rna1)[i - idx] = (*this)[i];
    }
    return rna1;
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
//-----------------------------------------------------------

RNA operator+(RNA &rna1, RNA &rna2){
    RNA result;
    result = rna1;
    for (size_t i = 1; i <= rna2.numb_of_nucl; i++){
        result.add_elem(rna2[i]);
    }
    return result;
}

bool RNA::operator==(const RNA &rna1) const{
    if (this->numb_of_nucl != rna1.numb_of_nucl) {
        return false;
    }
    if (!this->numb_of_nucl && !rna1.numb_of_nucl) {
        return true;
    }
    for (size_t i = 1; i < rna1.length_of_chain; i++){
        if (this->chain_of_nucl[i] != rna1.chain_of_nucl[i]) {
            return false;
        }
    }
    size_t idx = nucl_count_in_size_t * (length_of_chain - 1) + 1;
    for (idx; idx <= numb_of_nucl; idx++){
        if (this->GetNucl(idx) != rna1.GetNucl(idx)) {
            return false;
        }
    }
    return true;
}

bool RNA::operator!=(const RNA &rna1) const{//+
    return !((*this) == rna1);//((*this) == rnk) ? false : true;
}

RNA::reference RNA::operator[](size_t num){
    return reference(num, (*this));
}

Nucl RNA::operator[](size_t num) const{
    return GetNucl(num);
}

RNA RNA::operator!() const{//+
    RNA rna;
    rna.numb_of_nucl = this->numb_of_nucl;
    rna.length_of_chain = this->length_of_chain;
    rna.chain_of_nucl = new size_t[rna.length_of_chain];
    for (size_t i = 0; i < rna.length_of_chain; i++){
        rna.chain_of_nucl[i] = ~(this->chain_of_nucl[i]);
    }
    return rna;
}

RNA &RNA::operator=(const RNA &rna1){
    if (this == &rna1) {
        return *this;
    }
    if (this->length_of_chain != rna1.length_of_chain) {
        delete[] this->chain_of_nucl;
        this->chain_of_nucl = new size_t[rna1.length_of_chain];
        this->length_of_chain = rna1.length_of_chain;
    }

    this->numb_of_nucl = rna1.numb_of_nucl;
    for (size_t i = 0; i < this->length_of_chain; i++) {
        this->chain_of_nucl[i] = rna1.chain_of_nucl[i];
    }
    return (*this);
}

RNA::reference &RNA::reference::operator=(Nucl nucl){
    if (num > rna.numb_of_nucl){
        for (size_t i = 0; i < rna.numb_of_nucl - num; i++){
            rna.add_elem(nucl);
        }
    } else {
        size_t bitMask = 3;
        size_t idxArray = (size_t) ceil((float)num / (float)nucl_count_in_size_t) - 1;
        size_t shift = (nucl_count_in_size_t - (num - nucl_count_in_size_t * idxArray)) * 2;
        size_t n = (size_t) nucl << shift;
        rna.chain_of_nucl[idxArray] = (rna.chain_of_nucl[idxArray] & ~(bitMask << shift)) | n;
    }
    return (*this);
}

RNA::reference & RNA::reference::operator=(reference r){
    Nucl nucl = r.rna.GetNucl(r.num);
    rna[num] = nucl;
    return (*this);
}

RNA::reference::operator Nucl(){
    return rna.GetNucl(num);
}

//-----------------------------------------------------------

size_t RNA::length() const {
    return this->length_of_chain;
}

size_t RNA::num_of_nucls() const {
    return this->numb_of_nucl;
}

//void RNA::output() const {
//    int bits = sizeof(size_t) * 8 - 2;
//    size_t bitMask = 3 << bits;
//    size_t end_of_rna = 0;
//    for (size_t i = 0; i < length_of_chain; i++){
//        if (numb_of_nucl - i * nucl_count_in_size_t < nucl_count_in_size_t){
//            end_of_rna = numb_of_nucl - i * nucl_count_in_size_t;
//        }
//        else {
//            end_of_rna = nucl_count_in_size_t;
//        }
//        for (size_t j = 0; j < end_of_rna; j++){
//            int n = (int)((chain_of_nucl[i] & (bitMask >> (2 * j))) >> (bits - 2 * j));
//            switch (n){
//                case 0:
//                    std::cout << "A";
//                    break;
//                case 1:
//                    std::cout << "G";
//                    break;
//                case 2:
//                    std::cout << "C";
//                    break;
//                case 3:
//                    std::cout << "T";
//                    break;
//                default:
//                    break;
//            }
//        }
//        std::cout << " ";
//    }
//}



