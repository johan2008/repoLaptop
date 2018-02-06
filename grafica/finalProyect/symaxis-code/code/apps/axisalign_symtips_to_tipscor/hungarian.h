#ifndef _MY_HUNGARIAN_H
#define _MY_HUNGARIAN_H
#include <vector>
#include <iostream>
#include "matrix.h"
#include "munkres.h"
#include "VecDynamic.h"
using namespace std;
static void BipartiteHungarian(const vector< VecDynamic<double> >& data0, 
		const vector< VecDynamic<double> >& data1, vector<int>& ordered0, vector<int>& ordered1)
{
	Matrix<double> matrix(data0.size(), data1.size());
	ordered0.resize(data0.size());
	ordered1.resize(data1.size());
	for (int i=0; i<int(data0.size()); i++)
	{
		for (int j=0; j<int(data1.size()); j++)
		{
			const VecDynamic<double>& sample0 = data0[i];
			const VecDynamic<double>& sample1 = data1[j];
			VecDynamic<double> diff = sample0 - sample1;
			matrix(i, j) = diff.SumOfSquare();
		}
	}
	Munkres m;
	m.solve(matrix);
	for ( int row = 0 ; row < int(data0.size()) ; row++ ) {
		int rowcount = 0;
		ordered0[row] = row;
		ordered1[row] = -1;
		for ( int col = 0 ; col < int(data1.size()) ; col++  ) {
			if ( matrix(row,col) == 0 )
			{
				ordered1[row] = col;
				rowcount++;
			}
		}
		if ( rowcount != 1 )
			std::cerr << "Row " << row << " has " << rowcount << " columns that have been matched." << std::endl;
	}

	for ( int col = 0 ; col < int(data1.size()) ; col++ ) {
		int colcount = 0;
		for ( int row = 0 ; row < int(data0.size()) ; row++ ) {
			if ( matrix(row,col) == 0 )
				colcount++;
		}
		if ( colcount != 1 )
			std::cerr << "Column " << col << " has " << colcount << " rows that have been matched." << std::endl;
	}
}
#endif
