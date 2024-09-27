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
ThreeSizeInfo get3size(CGraph *gout, CGraph *gout_2) {

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
            
            for (EdgeIdx j = start; j < end; ++j) {
                const VertexIdx end1 = gout->nbors[j];
                
                // Tri4 : i < j < k (1, 1, 1)
                //        i->j : 1, i->k : 1, j->k : 1
                
                for (EdgeIdx k = j+1; k < end; ++k) {
                    const VertexIdx end2 = gout->nbors[k];
                    
                    EdgeIdx loc_111 = gout->getEdgeBinary(end1, end2);
                    if (loc_111 != -1) {
                        local_ret.tri4++;
                    }
                }

                // Tri2 : i < j < k 
                //    (1) i->j : 2, i->k : 1, j->k : 2
                //    (2) i->j : 1, i->k : 2, j->k : 2  

                for (EdgeIdx k = start_2; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];
                    EdgeIdx loc122 = end1 > end2 ? gout_2->getEdgeBinary(end2, end1) : gout_2->getEdgeBinary(end1, end2);
                    if (loc122 != -1) {
                        local_ret.tri2++;
                    }
                }                 
            }

            for (EdgeIdx j = start_2; j < end_2; ++j) {
                const VertexIdx end1 = gout_2->nbors[j];
                for (EdgeIdx k = j+1; k < end_2; ++k) {
                    const VertexIdx end2 = gout_2->nbors[k];

                    // Tri4 : i < j < k 
                    //    (3) i->j : 2, i->k : 2, j->k : 2

                    EdgeIdx loc_222 = gout_2->getEdgeBinary(end1, end2);
                    if (loc_222 != -1) {
                        local_ret.tri1++;
                    } else {
                        
                        // Tri2 : i < j < k 
                        //    (3) i->j : 2, i->k : 2, j->k : 1

                        EdgeIdx loc_221 = gout->getEdgeBinary(end1, end2);
                        if (loc_221 != -1) {
                            local_ret.tri2++;
                        }
                    }
                }
            }
            // #pragma omp critical
            // {
            //     if (current % 1000 == 0){
            //         printf("Node : %lld / %lld done... (node idx : %lld / out degree : %lld)\n", current, gout->nVertices, i, gout_2->degree(i));
            //     }
            //     current++;
            // }       
        }

        #pragma omp critical
        {
            ret.tri1 += local_ret.tri1;
            ret.tri2 += local_ret.tri2;
            ret.tri4 += local_ret.tri4;
        }
    }

    return ret;
}


void countThree(CGraph *cg, CDAG *dag, CGraph *cg_2, CDAG *dag_2, double (&mcounts)[6]){

    // 1. Count Three
    double w1 = 0, w2 = 0, t3 = 0; // # of 3-size d-Motifs : wedge 1, 2 & Triangle 1, 2, 3, 4 
    VertexIdx n = cg->nVertices;

    for (VertexIdx i = 0; i < n; i++){
        VertexIdx deg = cg->degree(i); // degree of i
        VertexIdx deg_2 = cg_2->degree(i); // degree2 of i

        w1 += deg_2 * (deg_2-1) / 2;
        w2 += deg * deg_2; 
        t3 += deg * (deg-1) / 2; 
    }

    mcounts[4] = w1;
    mcounts[5] = w2;
    mcounts[2] = t3;
    ThreeSizeInfo tricount = get3size(&(dag->outlist), &(dag_2->outlist));
    mcounts[0] = tricount.tri1;
    mcounts[1] = tricount.tri2;
    mcounts[3] = tricount.tri4;
}

void mEquation3(double (&mcounts)[6]){
    double t[6];

    t[0] = mcounts[0];
    t[1] = mcounts[1];
    t[2] = mcounts[2] - 3*mcounts[3];
    t[3] = mcounts[3];
    t[4] = mcounts[4] - 3 * t[0] - t[1];
    t[5] = mcounts[5] - 2 * t[1] - 2 * t[2];

    for (int m = 0; m < 36; ++m){
        mcounts[m] = t[m];
    }
}

void print3size(double (&mcounts)[6]){
    for (int i = 0; i < 6; ++i){
        printf("\"T%d\" : %.1f,\n", i+1, mcounts[i]);
    }
}
