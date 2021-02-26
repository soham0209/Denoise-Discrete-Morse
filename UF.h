#pragma once
#include <map>
#include <vector>
#include <iostream>
using std::vector;
class UF {
public:
	std::map<int, int> id;
	std::map<int, double> *f_val;
	int cnt;

	// Create an empty union find data structure with N isolated sets.
	UF(int N) {
		cnt = N;
		
	}
	~UF() {
	}
	void make(int p){
		id[p] = p;
	}
	// Return the id of component corresponding to object p.
	int root(int p) {
		while (p != id[p]){
			id[p] = id[id[p]];
			p = id[p];
		}
		return p;
	}
	int find(int p){
		return root(p);
	}	
	// Replace sets containing x and y with their union.
	void link(int i, int j) {
		//id[j] = i;
//		int rootp = root(i);
//        int rootq = root(j);
//        if (rootp == rootq) return;
//        double largestP = (*f_val)[large[rootp]];
//        double largestQ = (*f_val)[large[rootq]];
//
//		// make smaller root point to larger one
//		 if (sz[rootp] < sz[rootq]) {
//            id[rootp] = rootq;
//            sz[rootq] += sz[rootp];
//
//            if (largestP > largestQ)
//                large[rootq] = rootp;
//        } else {
//            id[rootq] = rootp;
//            sz[rootp] += sz[rootq];
//
//            if (largestQ > largestP)
//                large[rootp] = rootq;
//        }
    if ((*f_val)[i] > (*f_val)[j])
        id[j] = i;
    else if (std::abs((*f_val)[i] - (*f_val)[j]) < 1e-12) {
        if (i > j)
            id[j] = i;
        else
            id[i] = j;
    }
    else{
        id[i] = j;
    }
		cnt--;
	}
	// Are objects x and y in the same set?
	bool connected(int x, int y) {
		if (find(x) == -1 && find(y) == -1)
			return false;
		return find(x) == find(y);
	}
	// Return the number of disjoint sets.
	int count() {
		return cnt;
	}
	void merge_y_to_x(int x, int y) {
		id[y] = x;
		cnt--;
	}
};