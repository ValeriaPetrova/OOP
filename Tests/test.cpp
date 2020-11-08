#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "RNA.h"

TEST(TestRNAConstr, test1) {
    RNA* rna = new RNA();
    EXPECT_TRUE(rna != nullptr);
    EXPECT_TRUE(rna->length() == 0);
    EXPECT_TRUE(rna->num_of_nucls() == 0);
}

TEST(TestEquality, test2){
    RNA rna1;
    RNA rna2;
    EXPECT_EQ(true, rna1 == rna2);
}

TEST(TestNotEquality, test3){
    RNA rna1;
    RNA rna2;
    EXPECT_EQ(false, rna1 != rna2);
}

TEST(TestNegation, test4){
    RNA rna1;
    RNA rna2;
    rna2 = !rna1;
    EXPECT_EQ(true, rna1 == rna2);
}
