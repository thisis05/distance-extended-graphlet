#include "counter.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstdio> 
#include <filesystem>

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

    double mcounts[42];
    for (int m = 0; m < 42; ++m){
        mcounts[m] = 0;
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
            // double mcounts3[6];
            // printf("Count d-Motifs (3-size)\n");
            // auto start_time2 = high_resolution_clock::now();
            // countThree(&cg, &dag, &cg_2, &dag_2, mcounts3);
            // auto end_time2 = high_resolution_clock::now();
            // auto duration2 = duration_cast<milliseconds>(end_time2 - start_time2);
            // double seconds2 = duration2.count() / 1000.0; // Convert milliseconds to seconds
            // printf("Execution time for Counting Motifs (3-size): %.3f\n", seconds2);
            // mEquation3(mcounts3);

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

            //printf("Total Execution time for 3-size: %.3f\n", seconds1 + seconds2);
            printf("Total Execution time for 4-size: %.3f\n", seconds1 + seconds3);

            for (int m = 0; m < 42; ++m){
                if (m < 6)
                    continue;
                    //mcounts[m] += mcounts3[m];
                else
                    mcounts[m] += mcounts4[m-6];
            }
            file_count++;
        }
    }
    printf("===================================================\n");
    printf("Total Graphlet Average Count\n");
    // for (int i = 0; i < 6; ++i){
    //     printf("\"T%d\" : %.1f,\n", i+1, mcounts[i] / file_count);
    // }   

    for (int i = 6; i < 41; ++i){
        printf("\"Q%d\" : %.1f,\n", i-5, mcounts[i] / file_count);
    }
    printf("\"Q36\" : %.1f\n", mcounts[41] / file_count);
    return 0;
}