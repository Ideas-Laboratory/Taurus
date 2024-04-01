#ifndef _METRIC_H_
#define _METRIC_H_
#include <Node.h>
#include <Link.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <Edge.h>
#include <algorithm>

#define PI 3.141592653

double testSE(std::vector<Node> &nodes,std::vector<Link> &links,std::vector<std::vector<double> > &sp);
double testNP(vector<Node> &nodes,vector<Link> &links,vector<vector<double>> &sp,int neigh);
int judgeCross(Node x,Node y,Node z,Node h);
float multi(Node p1,Node p2,Node p0);
Node intersect(Node a,Node b,Node c,Node d);
double testCL(std::vector<Node> &nodes,std::vector<Link> &links);
double testMA(std::vector<Node> &nodes,std::vector<Link> &links);
double testCD(std::vector<Node> &nodes,std::vector<Link> &links,std::vector<int> &flag,int clusternum);
double testCE(std::vector<Node> &nodes, std::vector<Link> &links, std::vector<int> &flag, int clusternum);
#endif