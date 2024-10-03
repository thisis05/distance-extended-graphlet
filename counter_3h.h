#pragma once

#include "digraph.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)

struct ThreeSizeInfo {
    Count tri1; Count tri2; Count tri3; Count tri4;
    Count tri5; Count tri6; Count tri7; Count tri8;
    Count tri9; Count tri10; Count tri11; Count tri12;
    Count tri13;

    ThreeSizeInfo()
        : tri1(0), tri2(0), tri3(0), tri4(0), 
          tri5(0), tri6(0), tri7(0), tri8(0), 
          tri9(0), tri10(0), tri11(0), tri12(0), tri13(0){}

};

ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2, CGraph *gout_3);

void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, CGraph *cg_3, CDAG *dag_3, double (&mcounts)[13]);
void mEquation3(double (&mcounts)[13]);
void print3size(double (&mcounts)[13]);