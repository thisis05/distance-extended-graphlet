#include "counter_3h.h"
#include <algorithm>
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <stdexcept>
#include <cstring>
#ifdef _WIN32
    #include <windows.h>
    #include <psapi.h>
#endif
#include <iostream>
#include <chrono>

#define DEBUG_PRINT(fmt, ...) printf(fmt, __VA_ARGS__)

using namespace std;
using namespace std::chrono;

int num_threads = 6;
ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2, CGraph* gout_3) {

    omp_set_num_threads(num_threads);
    ThreeSizeInfo ret = {};
    VertexIdx current = 0;
    #pragma omp parallel
    {
        ThreeSizeInfo local_ret = {};
        
        #pragma omp for schedule(guided)
        for (VertexIdx i = 0; i < gout->nVertices; ++i) {
            const EdgeIdx start = gout->offsets[i];
            const EdgeIdx end = gout->offsets[i+1];
            const EdgeIdx start_2 = gout_2->offsets[i];
            const EdgeIdx end_2 = gout_2->offsets[i+1];
            const EdgeIdx start_3 = gout_3->offsets[i];
            const EdgeIdx end_3 = gout_3->offsets[i+1];
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                // 1 1 1
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.tri9++;
                    }
                }

                // 1 2 2
                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc122 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.tri7++;
                    }
                }       

                // 1 3 3
                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc133 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc133 != -1) {
                        local_ret.tri4++;
                    }
                }    

            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    // 2 2 3
                    EdgeIdx loc_223 = gout_3->getEdgeBinary(end1, end2);
                    if (loc_223 != -1){
                        local_ret.tri3++;
                    }
                    else{
                        // 2 2 2
                        EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                        if (loc_222 != -1) {
                            local_ret.tri6++;
                        } else {
                            // 2 2 1
                            EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                            if (loc_221 != -1) {
                                local_ret.tri7++;
                            }
                        }
                    }
                }
                for (EdgeIdx k = start_3; k < end_3; ++k) {
                    // 2 3 3
                    const VertexIdx end2 = gout_3->nbors[k];
                    EdgeIdx loc233 = end1 > end2 ? gout_3->getEdgeBinary(end2, end1) : gout_3->getEdgeBinary(end1, end2);
                    if (loc233 != -1) {
                        local_ret.tri2++;
                    }
                    else{
                        EdgeIdx loc232 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                        if (loc232 != -1) {
                            local_ret.tri3++;
                        }
                    }
                }
            }

            for (EdgeIdx j = start_3; j < end_3; ++j) {
                const VertexIdx end1 = gout_3->nbors[j];
                for (EdgeIdx k = j+1; k < end_3; ++k) {
                    const VertexIdx end2 = gout_3->nbors[k];
                    // 3 3 3
                    EdgeIdx loc_333 = gout_3->getEdgeBinary(end1, end2);
                    if (loc_333 != -1){
                        local_ret.tri1++;
                    }
                    else{
                        // 3 3 2
                        EdgeIdx loc_332 = gout_2->getEdgeBinary(end1, end2);
                        if (loc_332 != -1) {
                            local_ret.tri2++;
                        } else {
                            // 3 3 1
                            EdgeIdx loc_331 = gout->getEdgeBinary(end1, end2);
                            if (loc_331 != -1) {
                                local_ret.tri4++;
                            }
                        }
                    }
                }
            }
            #pragma omp critical
            {
                if (current % 100 == 0){
                    printf("Node : %lld / %lld done...\n", current, gout->nVertices);
                }
                current++;
            }       
        }

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri3 += local_ret.tri3;
            ret.tri4 += local_ret.tri4;
            ret.tri6 += local_ret.tri6;
            ret.tri7 += local_ret.tri7;
            ret.tri9 += local_ret.tri9;
        }
    }

    return ret;
}


void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, CGraph *cg_3, CDAG *dag_3, double (&mcounts)[13]){

    // 1. Count Three
    double t5 = 0, t8 = 0, t10 = 0, t11 = 0, t12 = 0, t13 = 0;
    VertexIdx n = cg->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = cg->degree(i); // degree of i
        VertexIdx deg_2 = cg_2->degree(i); // degree2 of i
        VertexIdx deg_3 = cg_3->degree(i); // degree2 of i

        t5 += deg * deg_2;
        t8 += deg * (deg-1) / 2; 
        t10 += deg_3 * (deg_3 - 1) / 2; 
        t11 += deg_2 * deg_3;
        t12 += deg_2 * (deg_2-1) / 2;
        t13 += deg * deg_3;
    }

    mcounts[4] = t5;
    mcounts[7] = t8;
    mcounts[9] = t10;
    mcounts[10] = t11;
    mcounts[11] = t12;
    mcounts[12] = t13;
    ThreeSizeInfo tricount = get3size(&(dag->outlist), &(dag_2->outlist), &(dag_3->outlist));
    mcounts[0] = tricount.tri1;
    mcounts[1] = tricount.tri2;
    mcounts[2] = tricount.tri3;
    mcounts[3] = tricount.tri4;
    mcounts[5] = tricount.tri6;
    mcounts[6] = tricount.tri7;
    mcounts[8] = tricount.tri9;
}

void mEquation3(double (&mcounts)[13]){
    double t[13];

    t[0] = mcounts[0];
    t[1] = mcounts[1];
    t[2] = mcounts[2];
    t[3] = mcounts[3];

    t[8] = mcounts[8];
    t[7] = mcounts[7] - 3 * t[8];
    t[6] = mcounts[6];

    t[4] = mcounts[4] - 2 * t[6] - 2 * t[7];
    t[5] = mcounts[5];
    
    t[9] = mcounts[9] - t[3] - t[1] - 3 * t[0];
    t[10] = mcounts[10] - 2 * t[1] - 2 * t[2] - t[4];
    t[11] = mcounts[11] - t[2] - 3 * t[5] - t[6];
    t[12] = mcounts[12] - 2 * t[3] - t[4];


    for (int m = 0; m < 13; ++m){
        mcounts[m] = t[m];
    }
}

void print3size(double (&mcounts)[13]){
    for (int i = 0; i < 13; ++i){
        printf("\"T%d\" : %.1f,\n", i+1, mcounts[i]);
    }
}
