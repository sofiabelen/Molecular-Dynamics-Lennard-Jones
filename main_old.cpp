#include <bits/stdc++.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <eggx.h>

using namespace std;

int Lx,Ly,Lz,win,n=125,m=5;
float kinetic,kinetic0,potential,potential0,maxv,minv,lx,ly,lz,mult=150.0;

struct vec
{
    float x,y,z;
};

vector<vec> pos(n),vel(n),acc(n),mv(27);

float randomRange(float a,float b)
{
    float random=(float(rand()))/float(RAND_MAX);
    return random*(b-a)+a;
}
void vecAdd(vec &a,vec b)
{
    a.x+=b.x;
    a.y+=b.y;
    a.z+=b.z;
}
void vecAdd(vec &a,vec b,float alpha)
{
    a.x+=b.x*alpha;
    a.y+=b.y*alpha;
    a.z+=b.z*alpha;
}
vec vecSub(vec a,vec b)
{
    vec c;
    c.x=a.x-b.x;
    c.y=a.y-b.y;
    c.z=a.z-b.z;
    return c;
}
vec vecScal(vec a,float s)
{
    a.x*=s;
    a.y*=s;
    a.z*=s;
    return a;
}
vec vecSet(float x,float y,float z)
{
    vec a;
    a.x=x;
    a.y=y;
    a.z=z;
    return a;
}
vec vecSet(float s)
{
    return vecSet(s,s,s);
}

void vecSet(vec &a,float x,float y,float z)
{
    a.x=x;
    a.y=y;
    a.z=z;
}

void vecSet(vec &a,float s)
{
    vecSet(a,s,s,s);
}
void vecRandNeg(vec &a)
{
    float tmp=rand();
    if(tmp>=RAND_MAX/2)
        a.x*=(-1);
    tmp=rand();
    if(tmp>=RAND_MAX/2)
        a.y*=(-1);
    tmp=rand();
    if(tmp>=RAND_MAX/2)
        a.z*=(-1);
}
float wrap(float a,float b)
{
    if(a<0)
        return a+float(int(int(-a)/int(b))+1)*b;
    else if(a>b)
        return int(a)%int(b-1)+float(int(a));
    else return a;
}
void vecWrap(vec &a)
{
    a.x=wrap(a.x,lx);
    a.y=wrap(a.y,ly);
    a.z=wrap(a.z,lz);
    //while(a.x>=lx)a.x-=lx;
    //while(a.y>=ly)a.y-=ly;
    //while(a.z>=lz)a.z-=lz;
    //while(a.x<0.0)a.x+=lx;
    //while(a.y<0.0)a.y+=ly;
    //while(a.z<0.0)a.z+=lz;
}
void vecWrapAll()
{
    for(int i=0;i<n;i++)
        vecWrap(pos[i]);
}
void vecMove(vec &a,vec v)
{
    vecAdd(a,v);
    vecWrap(a);
}
void redraw()
{
    gclr(win);
    for(int i=0;i<n;i++)
        fillcirc(win,pos[i].x*mult,pos[i].y*mult,5.0,5.0);
}
void vecUpdate()
{
    for(int i=0;i<n;i++)
        vecMove(pos[i],vel[i]);
    redraw();
}
float len(vec a)
{
    return sqrt(pow(a.x,2.0)+pow(a.y,2.0)+pow(a.z,2.0));
}
float vecDist(vec a,vec b)
{
    //vec r_ij=a;
    //vecAdd(r_ij,b,-1);
    //return len(r_ij);
    double r=INT_MAX;
    for(int j=0;j<27;j++)
    {
        vec b2=b;
        vecAdd(b2,mv[j]);
        vec r_ij=a;
        vecAdd(r_ij,b2,-1);
        r=min(double(r),double(len(r_ij)));
    }
    return r;
}
void energy()
{
    kinetic=0;
    potential=0;
    for(int i=0;i<n;i++)
    {
        minv=min(minv,sqrt(vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z));
        maxv=max(maxv,sqrt(vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z));
        kinetic+=vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z;
        for(int j=i+1;j<n;j++)
        {
            float d_ij=vecDist(pos[i],pos[j]);
            potential+=(pow(d_ij,-12)-pow(d_ij,-6));
            //printf("pot: %f\n",(pow(d_ij,-12)-pow(d_ij,-6)));
        }
       // printf("%.2f, ",vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z);
    }
    //printf("\n");
    //printf("%.3f %.3f\n",minv,maxv);
    //printf("E=%f, diff=%f\n",kinetic+potential,potential0+kinetic0-potential-kinetic);
    //printf("current value: %.3f starting value: %.3f difference: %.3f\n",kinetic,kinetic0,kinetic-kinetic0);
}
void init()
{
    float dx=lx/float(m);
    float dy=ly/float(m);
    float dz=lz/float(m);
    
    int u=0;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            for(int k=0;k<m;k++)
            {
                pos[u].x=randomRange((float(i)+0.4)*dx,(float(i)+0.6)*dx);
                pos[u].y=randomRange((float(j)+0.4)*dy,(float(j)+0.6)*dy);
                pos[u].z=randomRange((float(k)+0.4)*dz,(float(k)+0.6)*dz);

                //pos[u].x=(float(i)+0.5)*dx;
                //pos[u].y=(float(j)+0.5)*dy;
                //pos[u].z=(float(k)+0.5)*dz;

                vel[u].x=randomRange(0.00001,0.001);
                vel[u].y=randomRange(0.00001,0.001);
                vel[u].z=randomRange(0.00001,0.001);
                vecRandNeg(vel[u]);
                acc[u]=vecSet(0);
                u++;
            }
    energy();
    kinetic0=kinetic;
    potential0=potential;
    minv=float(INT_MAX);
    maxv=-float(INT_MAX);


    mv[0]=vecSet(-lx,0,0);
    mv[1]=vecSet(-lx,ly,0);
    mv[2]=vecSet(-lx,-ly,0);

    mv[3]=vecSet(-lx,0,lz);
    mv[4]=vecSet(-lx,ly,lz);
    mv[5]=vecSet(-lx,-ly,lz);
    
    mv[6]=vecSet(-lx,0,-lz);
    mv[7]=vecSet(-lx,ly,-lz);
    mv[8]=vecSet(-lx,-ly,-lz);

    mv[9]=vecSet(0,0,0);
    mv[10]=vecSet(0,ly,0);
    mv[11]=vecSet(0,-ly,0);

    mv[12]=vecSet(0,0,lz);
    mv[13]=vecSet(0,ly,lz);
    mv[14]=vecSet(0,-ly,lz);
    
    mv[15]=vecSet(0,0,-lz);
    mv[16]=vecSet(0,ly,-lz);
    mv[17]=vecSet(0,-ly,-lz);

    mv[18]=vecSet(lx,0,0);
    mv[19]=vecSet(lx,ly,0);
    mv[20]=vecSet(lx,-ly,0);

    mv[21]=vecSet(lx,0,lz);
    mv[22]=vecSet(lx,ly,lz);
    mv[23]=vecSet(lx,-ly,lz);
    
    mv[24]=vecSet(lx,0,-lz);
    mv[25]=vecSet(lx,ly,-lz);
    mv[26]=vecSet(lx,-ly,-lz);
}
void printVec(vec &a)
{
    printf("%.3f %.3f %.3f\n",a.x,a.y,a.z);
}
void printVec(vector<vec> &arr)
{
    for(int i=0;i<n;i++)
        printVec(arr[i]);
}
void forces()
{
    for(int i=0;i<n;i++)
    {
        vecAdd(vel[i],acc[i],1.0/2.0);
        vecAdd(pos[i],vel[i]);
        vecWrap(pos[i]);
    }
    for(int i=0;i<n;i++)
        acc[i]=vecSet(0);
    for(int i=0;i<n;i++)
    {
        for(int j=+1;j<n;j++)
        {
            //vec r_ij=add(pos[j],pos[i],-1.0);
            float d_ij=vecDist(pos[i],pos[j]);
            float f_ij=48*(pow(d_ij,-13)-0.5*(pow(d_ij,-7)));
            //printf("%f %f\n",d_ij,f_ij);
            float f_ijx=((pos[j].x-pos[i].x)/d_ij)*f_ij;
            float f_ijy=((pos[j].y-pos[i].y)/d_ij)*f_ij;
            float f_ijz=((pos[j].z-pos[i].z)/d_ij)*f_ij;
            acc[i].x+=f_ijx;
            acc[i].y+=f_ijx;
            acc[i].z+=f_ijx;
            acc[j].x-=f_ijx;
            acc[j].y-=f_ijx;
            acc[j].z-=f_ijx;
            
        }
    }
    for(int i=0;i<n;i++)
        printf("%f %f %f\n",acc[i].x,acc[i].y,acc[i].z);

    //printf("Accelerations\n");
    //printVec(acc);
    for(int i=0;i<n;i++)
        vecAdd(vel[i],acc[i],1.0/2.0);
}

int main()
{
    srand(time(NULL));
    ggetdisplayinfo(NULL,&Lx,&Ly);
    win=gopen(Lx,Ly);
    Lz=Lx;
    lx=float(Lx)/mult;
    ly=float(Ly)/mult;
    lz=float(Lz)/mult;
    gsetbgcolor(win,"BLACK");
    newcolor(win,"BLUE");
    //printf("density: %f\n",float(n)/(lx*ly*lz));

    init();
    //printf("mult:%.3f\n",mult);
    //printf("lx:%.3f ly:%.3f lz:%.3f\n",lx,ly,lz);
    //printf("Positions\n");
    //printVec(pos);
    //printf("Velocities\n");
    //printVec(vel);

    int i=0;
    while(i<1000)
    {
        forces();
        energy();
        redraw();
        msleep(1000);
        //printf("Positions\n");
        //printVec(pos);
        //printf("Velocities\n");
        //printVec(vel);
        //return 0;
        i++;
    }
    //gclose(win);
    return 0;
}
