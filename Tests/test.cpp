#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "RNA.h"

RNA newRNA(std::string rna){
    RNA rna1;
    for (size_t i = 0; i < rna.length(); i++){
        switch (rna[i]){
            case 'A':
                rna1.add_elem(A);
                break;
            case 'G':
                rna1.add_elem(G);
                break;
            case 'C':
                rna1.add_elem(C);
                break;
            case 'T':
                rna1.add_elem(T);
                break;
            default:
                break;
        }
    }
    return rna1;
}

TEST(TestRNA, test1) {
    RNA rna1(A, 10);
    RNA rna2 = newRNA("AAAAAAAAAA");
    EXPECT_TRUE(rna1 == rna2);
}

TEST(TestRNA, test2){
    RNA rna1;
    RNA rna2(rna1);
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestRNA, test3){
    RNA rna1 = newRNA("AAAACCCCGGGGTTTT");
    RNA rna2(rna1);
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestIsComplimentary, test4){
    RNA rna1 = newRNA("CCCCC");
    RNA rna2 = newRNA("GGGGG");
    EXPECT_EQ(true, rna1.isComplementary(rna2));
}

TEST(TestIsComplimentary, test5){
    RNA rna1 = newRNA("AAAAAAAA");
    RNA rna2 = newRNA("AAAAAAAA");
    EXPECT_EQ(false, rna1.isComplementary(rna2));
}

TEST(TestTrim,test6){
    RNA rna1 = newRNA("AGCTAGCTAGCT");
    RNA rna2 = newRNA("AGCT");
    rna1 = rna1.trim(4);
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestTrim, test31) {
    RNA rna1 = newRNA("AAAA");
    RNA rna2;
    rna1 = rna1.trim(0);
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestCardinality, test8){
    RNA rna(A, 10);
    rna[1] = rna[2] = rna[3] = T;
    EXPECT_EQ(rna.cardinality(T), 3);
}

TEST(TestCardinality, test9){
    RNA rna(A, 100);
    rna[1] = rna[2] = rna[3] = G;
    EXPECT_EQ(rna.cardinality(A), 100-4);
}

TEST(TestAdd, test10){
    RNA rna1 =newRNA("AAAAAAAA");
    RNA rna2;
    RNA rna3 = rna1 + rna2;
    RNA rna4 = rna2 + rna1;
    EXPECT_EQ(true, rna1 == rna3);
    EXPECT_EQ(true, rna1 == rna4);
}

TEST(TestAdd, test11){
    RNA rna1 = newRNA("AAAAAA");
    RNA rna2 = newRNA("CCCCCC");
    RNA rna3 = rna1 + rna2;
    RNA rna4 = newRNA("AAAAAACCCCCC");
    EXPECT_EQ(true, rna3 == rna4);
}

TEST(TestEqual, test26){
    RNA rna1;
    RNA rna2;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestEqual, test27){
    RNA rna1;
    RNA rna2 = newRNA("AAAA");
    EXPECT_EQ(false, rna1 == rna2);
}

TEST(TestEqual, test12){
    RNA rna1 = newRNA("AAAACCCCGGGGTTTT");
    RNA rna2 = newRNA("AAAACTTTGGGGCCCT");
    EXPECT_EQ(false, rna1 == rna2);
}

TEST(TestEqual, test13){
    RNA rna1 = newRNA("ACGTACGTACGT");
    RNA rna2 = newRNA("ACGTACGTACGT");
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestEqual, test14){
    RNA rna1 = newRNA("AAAA");
    RNA rna2 = newRNA("AAAAAAAA");
    EXPECT_EQ(false, rna1 == rna2);
}

TEST(TestNotEqual, test28){
    RNA rna1;
    RNA rna2;
    EXPECT_EQ(false, rna1 != rna2);
}

TEST(TestNotEqual, test29){
    RNA rna1;
    RNA rna2 = newRNA("AAAA");
    EXPECT_EQ(true, rna1 != rna2);
}

TEST(TestNotEqual, test15){
    RNA rna1 = newRNA("ACGTACGTACGT");
    RNA rna2 = newRNA("ACGTACGTACGT");
    EXPECT_EQ(false, rna1 != rna2);
}

TEST(TestNotEqual, test16){
    RNA rna1 = newRNA("AAAACCCCGGGGTTTT");
    RNA rna2 = newRNA("AAAACTTTGGGGCCCT");
    EXPECT_EQ(true, rna1 != rna2);
}

TEST(TestNotEqual, test17){
    RNA rna1 = newRNA("AAAA");
    RNA rna2 = newRNA("AAAAAAAA");
    EXPECT_EQ(true, rna1 != rna2);
}

TEST(TestNot, test18){
    RNA rna1;
    RNA rna2 = newRNA("GCTG");
    rna1 = !rna2;
    EXPECT_EQ(false, rna1 == rna2);
}

TEST(TestNot, test19){
    RNA rna1 = newRNA("CGAC");
    RNA rna2 = newRNA("GCTG");
    rna1 = !rna1;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestNot, test20){
    RNA rna1 = newRNA("CGAC");
    RNA rna2 = newRNA("CGAC");
    rna2 = !rna2;
    EXPECT_EQ(false, rna1 == rna2);
}

TEST(TestNot, test30){
    RNA rna1;
    RNA rna2;
    rna2 = !rna1;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestAllocation, test21){
    RNA rna1 = newRNA("AAAACTTTGGGGCCCT");
    RNA rna2;
    rna2 = rna1;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestAllocation, test22){
    RNA rna1;
    RNA rna2 = newRNA("AAAACCCCGGGGTTTT");
    rna1 = rna2;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestAllocation, test23){
    RNA rna1 = newRNA("ACGTACGTACGT");
    RNA rna2 = newRNA("AAAACTTTGGGGCCCT");
    rna2 = rna1;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestAllocation, test29){
    RNA rna1;
    RNA rna2;
    rna2 = rna1;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(LargeTest, test24){
    RNA rna1(A, 1000000);
    RNA rna2;
    for (size_t i = 1; i <= 1000000; i++){
        rna2[i] = T;
    }

    EXPECT_EQ(true, rna1 != rna2);
}

TEST(TestLength, test25){
    RNA rna(A, 10000);
    size_t allocLength = rna.length();
    EXPECT_EQ(allocLength, 10000 * 2 / 8 / sizeof(size_t) + 1);
}
