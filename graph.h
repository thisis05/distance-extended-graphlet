#pragma once

#include <cstdint>
#include <map>
#include <vector>
#include <unordered_map>
#include <string>

using VertexIdx = uint32_t;
using EdgeIdx = uint64_t;
using Count = uint64_t;

struct Graph
{
  VertexIdx  nVertices; //number of vertices in the graph
  EdgeIdx    nEdges;    //number of edges in this list
  VertexIdx *srcs;      //array of source vertices
  VertexIdx *dsts;      //array of destination vertices
  
  Graph copy() const;
};

Graph newGraph(VertexIdx nVertices, EdgeIdx nEdges);

struct CGraph
{
  VertexIdx  nVertices;   //number of vertices in the graph
  EdgeIdx    nEdges;      //number of edges in the graph
  EdgeIdx   *offsets;     //vertex v's edges are [offset[v], offset[v + 1])
  VertexIdx *nbors;       //incoming or outgoing neighbors
  VertexIdx *mapping;
  VertexIdx *inverse;
  Count maxDegree;
  //std::unordered_map<long long, std::unordered_map<long long, long long>> adj_matrix;
  
  CGraph renameByDegreeOrder() const;
  CGraph reMapping(VertexIdx *mapping, VertexIdx *inverse) const;
  void sortById() const;
  CGraph getE2() const;
  CGraph getE3(CGraph G1) const;
  CGraph getE4(CGraph G1, CGraph G2) const;
  EdgeIdx getEdgeBinary(VertexIdx v1, VertexIdx v2) const;
  EdgeIdx degree(VertexIdx v) const {return offsets[v + 1] - offsets[v]; }
  bool isEdgeBinary(VertexIdx v1, VertexIdx v2) const;
};

CGraph makeCSR(Graph& graph);
CGraph newCGraph(VertexIdx nVertices, EdgeIdx nEdges);

struct Pair {
    VertexIdx first;
    VertexIdx second;
};

bool compareDegree(Pair firstPair, Pair nextPair);
VertexIdx binarySearch(EdgeIdx* array, VertexIdx end, EdgeIdx val);
void read_mtx(const std::string& filename, Graph& graph);

