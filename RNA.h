#ifndef RNA_RNA_H
#define RNA_RNA_H

#include <iostream>
#include <cstdint>
#include <cmath>

enum Nucl{
    A, //0 == 00
    G, //1 == 01
    C, //2 == 10
    T  //3 == 11
};

const size_t nucl_count_in_size_t = sizeof(size_t) * 8 / 2;

class RNA{
private:
    size_t *chain_of_nucl = nullptr; // цепь нуклеатидов
    size_t length_of_chain = 0; // количество size_t
    size_t numb_of_nucl = 0; // количество нуклеатидов
public:
    //____________________
    class reference{
    private:
        size_t num;
        RNA &rna;
    public:
        reference(size_t, RNA&);

        reference &operator=(Nucl);
        reference &operator=(reference);
        operator Nucl();

        ~reference();
    };
    //______________________
    RNA();
    RNA(Nucl, size_t);
    RNA(const RNA&);

    virtual ~RNA();// деструктор

    Nucl GetNucl(size_t) const;
    void add_elem(Nucl);
    Nucl complementary(Nucl) const;
    bool isComplementary(RNA&) const;
    RNA& trim(size_t);
    RNA& split(size_t);
    size_t cardinality(Nucl);

    friend RNA operator+ (RNA&, RNA&);
    bool operator== (const RNA&) const;
    bool operator!= (const RNA&) const;

    reference operator[] (size_t);
    Nucl operator[](size_t) const;
    RNA operator! () const;
    RNA& operator= (const RNA&);

//    void output() const;
    size_t length() const;
    size_t num_of_nucls() const;
};

#endif //RNA_RNA_H
