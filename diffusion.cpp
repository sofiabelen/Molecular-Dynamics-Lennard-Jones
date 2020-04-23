#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

ifstream positions;
ofstream diff;
int n, ntime;
double dt;
vector<vector<vector<double> > > pos;
vector<double> msd;

int main()
{
    positions.open("positions");
    diff.open("diff");
    diff<<"t x"<<endl;
    positions>>n>>ntime>>dt;

    pos.resize(ntime,vector<vector<double> >(n,vector<double>(3)));
    
    for(int t=0;t<ntime;t++)
        for(int i=0;i<n;i++)
            for(int j=0;j<3;j++)
                positions>>pos[t][i][j];

    int nmin = int((ntime-1)/2);
    msd.resize(nmin,0);

    //Iterate over all delta T
    for(int k=1;k<=nmin;k++)
    {
        //Iterate over x,y,z
        for(int d=0;d<3;d++)
        {
            //Iterate over all molecules
            for(int i=0;i<n;i++)
            {
                //Iterate over all time origins
                for(int t0=0;t0<nmin;t0++)
                {
                    msd[k-1] += (pos[t0+k][i][d]-pos[t0][i][d])*
                                (pos[t0+k][i][d]-pos[t0][i][d]);
                }
            }
            msd[k-1] /= double(n*nmin*k)*3.0;
        }
        diff<<k*dt<<" "<<msd[k-1]<<endl;
    }
    return 0;
}
