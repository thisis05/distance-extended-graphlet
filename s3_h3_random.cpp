#include "counter_3h.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstdio> 

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

//#define PRINT_CSR

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string directory_path = argv[1];
    int file_count = 0;

    double mcounts[13];
    for (int m = 0; m < 13; ++m){
        mcounts[m] = 13;
    };

    for (const auto& entry : fs::directory_iterator(directory_path)) {
        printf("Random Graph %d \n", file_count+1);
        if (entry.is_regular_file() && entry.path().extension() == ".mtx") {
            string filename = entry.path().string();

            Graph g;
            printf("Read Graph\n");
            read_mtx(filename, g);

            printf("Loaded graph\n");
            CGraph pre_cg = makeCSR(g);
            pre_cg.sortById();
            printf("Converted to CSR\n");

            auto start_time1 = high_resolution_clock::now();

            printf("Get Edge set with distance 2\n");
            CGraph pre_cg_2 = pre_cg.getE2();
            pre_cg_2.sortById();

            auto end_time1 = high_resolution_clock::now();
            auto duration1 = duration_cast<milliseconds>(end_time1 - start_time1);
            double seconds1 = duration1.count() / 1000.0; // Convert milliseconds to seconds
            printf("Execution time for get E2:  %.3f\n", seconds1);

            auto start_time1_1 = high_resolution_clock::now();

            printf("Get Edge set with distance 3\n");
            CGraph pre_cg_3 = pre_cg_2.getE3(pre_cg);
            
            auto end_time1_1 = high_resolution_clock::now();
            auto duration1_1 = duration_cast<milliseconds>(end_time1_1 - start_time1_1);
            double seconds1_1 = duration1_1.count() / 1000.0; // Convert milliseconds to seconds
            printf("Execution time for get E3:  %.3f\n", seconds1_1);
    
            printf("Creating DAG\n");
            CGraph cg_3 = pre_cg_3.renameByDegreeOrder();
            cg_3.sortById();
            CGraph cg = pre_cg.reMapping(cg_3.mapping, cg_3.inverse);
            CGraph cg_2 = pre_cg_2.reMapping(cg_3.mapping, cg_3.inverse);
            cg.sortById();
            cg_2.sortById();

            CDAG dag_3 = degreeOrdered(&cg_3);
            (dag_3.outlist).sortById();
            (dag_3.inlist).sortById();

            CDAG dag_2 = degreeOrdered2(&cg_2, &cg_3);
            (dag_2.outlist).sortById();
            (dag_2.inlist).sortById();

            CDAG dag = degreeOrdered2(&cg, &cg_3);
            (dag.outlist).sortById();
            (dag.inlist).sortById();

            
            //1. Count 3-size d-Motifs
            double mcounts3[13];
            printf("Count d-Motifs (3-size)\n");
            auto start_time2 = high_resolution_clock::now();
            countThree(&cg, &dag, &cg_2, &dag_2, &cg_3, &dag_3, mcounts3);
            auto end_time2 = high_resolution_clock::now();
            auto duration2 = duration_cast<milliseconds>(end_time2 - start_time2);
            double seconds2 = duration2.count() / 1000.0; // Convert milliseconds to seconds
            printf("Done\n");
            printf("Execution time for Counting Motifs (3-size): %.3f\n", seconds2);
            mEquation3(mcounts3);
            print3size(mcounts3);

            printf("Total Execution time for 3-size: %.3f\n", seconds1 + seconds2);

            for (int m = 0; m < 13; ++m){
                mcounts[m] += mcounts3[m];
            }
            file_count++;
    
        }
    }
    printf("===================================================\n");
    printf("Total Graphlet Average Count\n");
    for (int i = 0; i < 12; ++i){
        printf("\"T%d\" : %.1f,\n", i+1, mcounts[i] / file_count);
    }   
    printf("\"T13\" : %.1f\n", mcounts[12] / file_count);
    return 0;
}