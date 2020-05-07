#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

ofstream diff;
ifstream tracker, positions;
int n, ntime, exp_count;
double dt;
vector<vector<vector<double> > > pos;
vector<double> msd,l;

int main()
{
    exp_count=-1;
    tracker.open("Data/counter");
    tracker>>exp_count;
    tracker.close();

    exp_count--;
    positions.open("Data/positions"+to_string(exp_count));
    diff.open("Data/diff"+to_string(exp_count));
    diff<<"t x"<<endl;
    l.resize(3);
    positions>>n>>ntime>>dt;

    for(int i=0;i<3;i++)
        positions>>l[i];

    pos.resize(ntime,vector<vector<double> >(n,vector<double>(3)));
    
    for(int t=0;t<ntime;t++)
        for(int i=0;i<n;i++)
            for(int j=0;j<3;j++)
                positions>>pos[t][i][j];

    vector<double> add(3,0);
    //for(int t=1;t<ntime;t++)
    //    for(int i=0;i<n;i++)
    //        for(int j=0;j<3;j++)
    //        {
    //            if(pos[t-1][i][j] - pos[t][i][j] <= l[j]/2.0)
    //                add[j] -= l[j];
    //            if(pos[t-1][i][j] - pos[t][i][j] >= l[j]/2.0)
    //                add[j] += l[j];
    //            pos[t][i][j] += add[j];
    //        }

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
            msd[k-1] /= double(n*nmin)*3.0;
        }
        diff<<k*dt<<" "<<msd[k-1]<<endl;
    }
    return 0;
}
