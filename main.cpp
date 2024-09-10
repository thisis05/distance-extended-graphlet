#include "counter.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstdio> 

using namespace std;
using namespace std::chrono;

//#define PRINT_CSR

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    string filename = argv[1];
    Graph g;
    printf("Read Graph\n");
    read_mtx(filename, g);

    printf("Loaded graph\n");
    CGraph pre_cg = makeCSR(g);
    printf("Converted to CSR\n");

    auto start_time1 = high_resolution_clock::now();

    printf("Get Edge set with distance 2\n");
    CGraph pre_cg_2 = pre_cg.getE2();
    CGraph cg_2 = pre_cg_2.renameByDegreeOrder();
    cg_2.sortById();

    auto end_time1 = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(end_time1 - start_time1);
    double seconds1 = duration1.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for get E2:  %.3f\n", seconds1);

    printf("Creating DAG\n");
    CGraph cg = pre_cg.reMapping(cg_2.mapping, cg_2.inverse);
    cg.sortById();

    CDAG dag_2 = degreeOrdered(&cg_2);
    (dag_2.outlist).sortById();
    (dag_2.inlist).sortById();

    CDAG dag = degreeOrdered2(&cg, &cg_2);
    (dag.outlist).sortById();
    (dag.inlist).sortById();

    
    // 1. Count 3-size d-Motifs
    double mcounts3[6];
    printf("Count d-Motifs (3-size)\n");
    auto start_time2 = high_resolution_clock::now();
    countThree(&cg, &dag, &cg_2, &dag_2, mcounts3);
    auto end_time2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(end_time2 - start_time2);
    double seconds2 = duration2.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting Motifs (3-size): %.3f\n", seconds2);
    mEquation3(mcounts3);

    // 2. Count 4-size d-Motifs
    double mcounts4[36];
    printf("Count d-Motifs (4-size)\n");
    auto start_time3 = high_resolution_clock::now();
    countFour(&dag, &dag_2, mcounts4);
    auto end_time3 = high_resolution_clock::now();
    auto duration3 = duration_cast<milliseconds>(end_time3 - start_time3);
    double seconds3 = duration3.count() / 1000.0; // Convert milliseconds to seconds
    printf("Execution time for Counting Motifs (4-size): %.3f\n", seconds3);
    mEquation4(mcounts4);

    printf("Total Execution time for 3-size: %.3f\n", seconds1 + seconds2);
    printf("Total Execution time for 4-size: %.3f\n", seconds1 + seconds3);

    printf("# of Edge (1): %lld", cg.nEdges / 2);
    printf("\n# of Edge (2): %lld", cg_2.nEdges / 2);
    printf("\nMax Degree (1): %lld", cg.maxDegree);
    printf("\nMax Degree (2): %lld\n", cg_2.maxDegree);

    
    //FILE* file = fopen("./degree/temp.txt", "w");
    VertexIdx degree = 0, degree2 = 0, degree3 = 0, degree4 = 0; 
    Count degree_sum = 0, degree2_sum = 0, degree3_sum = 0, degree4_sum = 0;

    //fprintf(file, "outdeg1 indeg1 totaldeg1 outdeg2 indeg2 totaldeg2 totaldeg\n");
    for (VertexIdx i = 0; i < dag.outlist.nVertices; ++i) {
        VertexIdx off = dag.outlist.offsets[i+1] - dag.outlist.offsets[i];
        degree_sum += off;
        if (off > degree){
            degree = off;
        }

        VertexIdx off2 = dag.inlist.offsets[i+1] - dag.inlist.offsets[i];
        degree2_sum += off2;
        if (off2 > degree2){
            degree2 = off2;
        }

        VertexIdx off3 = dag_2.outlist.offsets[i+1] - dag_2.outlist.offsets[i];
        degree3_sum += off3;
        if (off3 > degree3){
            degree3 = off3;
        }

        VertexIdx off4 = dag_2.inlist.offsets[i+1] - dag_2.inlist.offsets[i];
        degree4_sum += off4;
        if (off4 > degree4){
            degree4 = off4;
        }

        VertexIdx total_off1 = off + off2;
        VertexIdx total_off2 = off3 + off4;
        VertexIdx total_off = total_off1 + total_off2;
        //fprintf(file, "%lld %lld %lld %lld %lld %lld %lld\n", off, off2, total_off1, off3, off4, total_off2, total_off);
    }
    //fclose(file);
    printf("DAG (Out) ______________________________________________________\n");
    
    printf("Average Degree : %.lld\n", degree_sum / dag.outlist.nVertices);
    printf("Max Degreee : %lld\n", degree);

    printf("DAG (in) ______________________________________________________\n");
   
    printf("Average Degree : %.lld\n", degree2_sum / dag.inlist.nVertices);
    printf("Max Degree : %lld\n", degree2);

    printf("DAG 2 (Out) ______________________________________________________\n");

    printf("Average Degree : %.lld\n", degree3_sum / dag_2.outlist.nVertices);
    printf("Max Degree : %lld\n", degree3);

    printf("DAG 2 (in) ______________________________________________________\n");

    printf("Average Degree : %.lld\n", degree4_sum / dag_2.inlist.nVertices);
    printf("Max Degree : %lld\n", degree4);


    #ifdef PRINT_CSR

    printf("CGraph ______________________________________________________\n");
    printf("Offsets ");
    for (VertexIdx i = 0; i <= cg.nVertices; ++i) {
        printf("%lld ", cg.offsets[i]);
    }
    printf("\n");


    printf("Destinations : ");
    for (EdgeIdx i = 0; i < cg.nEdges; ++i) {
        printf("%lld ", cg.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < cg.nVertices; ++i) {
        printf("%lld ", cg.degree(i));
    }
    printf("\n");
    
    printf("Offsets ");
    for (VertexIdx i = 0; i <= cg_2.nVertices; ++i) {
        printf("%lld ", cg_2.offsets[i]);
    }
    printf("\n");

    printf("Destinations : ");
    for (EdgeIdx i = 0; i < cg_2.nEdges; ++i) {
        printf("%lld ", cg_2.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < cg_2.nVertices; ++i) {
        printf("%lld ", cg_2.degree(i));
    }
    printf("\n");

    printf("DAG (Out) ______________________________________________________\n");
    printf("Offsets ");
    for (VertexIdx i = 0; i <= dag.outlist.nVertices; ++i) {
        printf("%lld ", dag.outlist.offsets[i]);
    }
    printf("\n");


    printf("Destinations : ");
    for (EdgeIdx i = 0; i < dag.outlist.nEdges; ++i) {
        printf("%lld ", dag.outlist.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < dag.outlist.nVertices; ++i) {
        printf("%lld ", dag.outlist.degree(i));
    }
    printf("\n");

    printf("Offsets ");
    for (VertexIdx i = 0; i <= dag_2.outlist.nVertices; ++i) {
        printf("%lld ", dag_2.outlist.offsets[i]);
    }
    printf("\n");

    printf("Destinations : ");
    for (EdgeIdx i = 0; i < dag_2.outlist.nEdges; ++i) {
        printf("%lld ", dag_2.outlist.nbors[i]);
    }
    printf("\n");

    printf("Degree: ");
    for (VertexIdx i = 0; i < dag_2.outlist.nVertices; ++i) {
        printf("%lld ", dag_2.outlist.degree(i));
    }
    printf("\n");
    #endif
    
    return 0;
}