#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <set>
#include <string>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

const int threshold = 2;

tuple<vector<vector<int>>, vector<vector<int>>> readGraph(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Failed to open file: " << filename << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        if (line[0] != '%') break;
    }

    int n, m;
    istringstream iss(line);
    iss >> n >> n >> m;

    printf("Read File ... \n");
    vector<vector<int>> adjMatrix(n, vector<int>(n,0));
    vector<vector<int>> adjVector(n);

    int u, v;
    while (file >> u >> v) {
        u--, v--;
        adjMatrix[u][v] = 1;
        adjMatrix[v][u] = 1;  
        adjVector[u].push_back(v);
        adjVector[v].push_back(u);
    }
    file.close();


    printf("Calculate Distance ... \n");
    vector<vector<int>> dist = adjMatrix;

    for (int start_node = 0; start_node < n; start_node++) {
        vector<int> second = adjVector[start_node];
        for (const auto& second_node : second) {
            vector<int> final = adjVector[second_node];
            for (const auto& final_node : final) {
                if (start_node != final_node && dist[start_node][final_node] != 1){
                    dist[start_node][final_node] = 2;
                    adjMatrix[start_node][final_node] = 1;
                }
            }
        }
    }

    int count_e3 = 0;
    for (int start_node = 0; start_node < n; start_node++) {
        vector<int> second = adjVector[start_node];
        for (const auto& second_node : second) {
            vector<int> third = adjVector[second_node];
            for (const auto& third_node : third) {
                if (start_node != third_node){
                    vector<int> final = adjVector[third_node];
                    for (const auto& final_node : final) {
                        if (start_node != final_node && second_node != final_node
                            && dist[start_node][final_node] != 1 && dist[start_node][final_node] != 2){
                            if (dist[start_node][final_node] != 3){
                                dist[start_node][final_node] = 3;
                                adjMatrix[start_node][final_node] = 1;
                                count_e3++;
                            }
                        }
                    }
                }
            }
        }
    }
    printf("E3 : %d\n", count_e3 / 2);

    return make_tuple(dist, adjMatrix);
}


void GenerateSubgraphsFromEdge(int u, int v, const vector<unordered_set<int>>& adjList, const vector<vector<int>>& adjMatrix, set<int>& extension) {
    // u의 이웃 노드를 extension에 추가
    for (int neighbor : adjList[u]) {
        if (neighbor != v) {
            extension.insert(neighbor);
        }
    }

    // v의 이웃 노드를 extension에 추가
    for (int neighbor : adjList[v]) {
        if (neighbor != u) {
            extension.insert(neighbor);
        }
    }
}

map<string, long long> countMotifs(vector<vector<int>> dist, const vector<vector<int>>& adjMatrix) {

    int n = adjMatrix.size();
    vector<pair<int, int>> edges;

    vector<unordered_set<int>> adjList(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (adjMatrix[i][j] > 0 && i != j) { 
                adjList[i].insert(j);
                if (i < j)
                    edges.emplace_back(i, j);
            }
        }
    }

    map <string, long long> motifCounts = {
        {"T1", 0}, {"T2", 0}, {"T3", 0}, {"T4", 0}, {"T5", 0}, {"T6", 0}, {"T7", 0}, 
        {"T8", 0}, {"T9", 0}, {"T10", 0}, {"T11", 0}, {"T12", 0}, {"T13", 0}, 
    };

    set<int> T1 = {6, 6, 6}; set<int> T2 = {5, 5, 6}; set<int> T3 = {4, 5, 5}; set<int> T4 = {4, 4, 6}; 
    set<int> T5 = {3, 4, 5}; set<int> T6 = {4, 4, 4}; set<int> T7 = {3, 3, 4}; set<int> T8 = {2, 3, 3};
    set<int> T9 = {2, 2, 2}; set<int> T10 = {3, 3, 6}; set<int> T11 = {2, 3, 5}; set<int> T12 = {2, 2, 4};
    set<int> T13 = {1, 3, 4};

    int max_num_workers = omp_get_max_threads()-1;
    omp_set_num_threads(max_num_workers);
    printf("Get Max Thread : %d\n", max_num_workers);    
    vector<map<string, long long>> thread_motifCounts(max_num_workers, motifCounts);
    int total_new_edges = edges.size();
    int progress = 0;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        #pragma omp for schedule(dynamic)
        // 1st Loop: Iterate over each key in the map
        for (int i = 0; i < edges.size(); ++i) {

            int a = edges[i].first;
            int b = edges[i].second;
            set<int> extension;

            GenerateSubgraphsFromEdge(a, b, adjList, adjMatrix, extension);

            for (int c : extension) {
                int structure = (adjMatrix[a][b] + adjMatrix[a][c] + adjMatrix[b][c]);
                int dist_ab = dist[a][b];
                int dist_ac = dist[a][c];
                int dist_bc = dist[b][c];
                int total_degree = 2 * (dist_ab + dist_ac + dist_bc);
                set<int> degree;
                degree.insert(dist_ab + dist_ac); degree.insert(dist_ab + dist_bc); degree.insert(dist_ac + dist_bc);

                if (structure == 2){
                    if (total_degree == 8){
                        if(degree == T13)
                            thread_motifCounts[thread_id]["T13"]++;
                        else
                            thread_motifCounts[thread_id]["T12"]++;
                    }
                    else if(total_degree == 10){
                        thread_motifCounts[thread_id]["T11"]++;
                    }
                    else if(total_degree == 12){
                        thread_motifCounts[thread_id]["T10"]++;
                    }
                    else{
                        printf("remain star type : %d\n", total_degree);
                    }
                }
                else if (structure == 3){
                    if (total_degree == 6)
                        thread_motifCounts[thread_id]["T9"]++;
                    else if (total_degree == 8)
                        thread_motifCounts[thread_id]["T8"]++;
                    else if (total_degree == 10)
                        thread_motifCounts[thread_id]["T7"]++;
                    else if (total_degree == 12){
                        if(degree == T6)
                            thread_motifCounts[thread_id]["T6"]++;
                        else
                            thread_motifCounts[thread_id]["T5"]++;
                    }
                    else if (total_degree == 14){
                        if(degree == T4)
                            thread_motifCounts[thread_id]["T4"]++;
                        else
                            thread_motifCounts[thread_id]["T3"]++;
                    }
                    else if (total_degree == 16)
                        thread_motifCounts[thread_id]["T2"]++;
                    else if (total_degree == 18)
                        thread_motifCounts[thread_id]["T1"]++;
                    else{
                        printf("remain triangle type : %d\n", total_degree);
                    }
                }
                else{
                    printf("remain structure : %d\n", structure);
                }
            }
            #pragma omp atomic
            progress++;
            if (progress % 100 == 0) {
                #pragma omp critical
                {
                    printf("Progress: %d / %d\n", progress, total_new_edges);
                }
            }
        }
    }
    for (int i = 0; i < max_num_workers; i++) {
        for (const auto& motif : thread_motifCounts[i]) {
            motifCounts[motif.first] += motif.second;
        }
    }
    return motifCounts;
}


void printAdjMatrix(vector<vector<int>> adjMatrix) {
    int n = adjMatrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

void printDistMatrix(vector<vector<int>> dist) {
    int n = dist.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << dist[i][j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    auto start_time = high_resolution_clock::now();

    string filename = argv[1];
    tuple<vector<vector<int>>, vector<vector<int>>> result = readGraph(filename);
    vector<vector<int>> dist = get<0>(result);
    vector<vector<int>> adjMatrix = get<1>(result);

    // cout << "Adjacency Matrix:" << endl;
    // printAdjMatrix(adjMatrix);
    // cout << "Dist Matrix:" << endl;
    // printDistMatrix(dist);
    map <string, long long> results = countMotifs(dist, adjMatrix);

    for (auto& motif : results) {
        if (motif.first[0] == 'T') {
            // Extract the number part after "T"
            int num = stoi(motif.first.substr(1));  // "T" 이후의 숫자 부분 추출
            if (num >= 1 && num <= 9) {
                motif.second /= 3;
            } else if (num >= 10 && num <= 13) {
                motif.second /= 2;
            }
        }
    }
    

    cout << "Motif Counts:" << endl;
    for (const auto& motif : results) {
        cout << "\"" << motif.first << "\"" << " : " << motif.second << "," << endl;
    }
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end_time - start_time);
    cout << "Execution time: " << duration.count() << " seconds" << endl;

    return 0;
}