//
// Created by VALERIA on 02.11.2020.
//
#ifndef UNTITLED2_RNA_H
#define UNTITLED2_RNA_H

#include <cstdio>
#include <cstdlib>

enum Nucl{
    A, //0 == 00
    G, //1 == 01
    C, //2 == 10
    T  //3 == 11
};

class RNA{
private:
    size_t *chain_of_nucl; // цепь нуклеатидов
    size_t length_of_chain; // количество size_t
    size_t numb_of_nucl; // количество нуклеатидов
    size_t nucl_count_in_size_t = sizeof(size_t) * 8 / 2;
public:
    //____________________
    class reference{
    private:
        size_t num;
        RNA *rna;
    public:
        reference(size_t, RNA*);

        reference &operator=(Nucl n);
        operator Nucl() const;
        ~reference();
    };
    //______________________
    RNA() : chain_of_nucl(nullptr), numb_of_nucl(0), length_of_chain(0) { };
    RNA(Nucl N, size_t length);
    explicit RNA(size_t);
    RNA(const RNA&);

    virtual ~RNA();// деструктор

    Nucl GetNucl(size_t);
    void add_elem(Nucl);
    Nucl complementary(Nucl);
    bool isComplementary(Nucl, Nucl);
    void trim(size_t);
    RNA split(size_t index);
    size_t cardinality(Nucl);

    friend RNA operator+ (RNA&, RNA&);
    bool operator== (const RNA&) const;
    bool operator!= (const RNA&) const;
    RNA operator! ();
    RNA& operator= (RNA const &);
    reference operator[] (size_t num);

    void output() const;
};


#endif //UNTITLED2_CLASSNAME_H
