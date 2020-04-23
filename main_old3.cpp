#include <bits/stdc++.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <eggx.h>

using namespace std;

int Lx,Ly,Lz,win,m,counter,border,border2,N,n,win2,l_hist,d,h_frames,wait_time;
double kinetic,kinetic0,potential,potential0,lx,ly,lz,mult,dt,density,pressure,vir_sum,v_max,temperature;
vector<int> x_vel;
string bglight,bgdark,accent,primary,divider,darkprimary,light;

struct vec
{
    double x,y,z;
};

vector<vec> pos(n),vel(n),acc(n),mv(27);

double randomRange(double a,double b)
{
    double random=(double(rand()))/double(RAND_MAX);
    return random*(b-a)+a;
}
void vecAdd(vec &a,vec b)
{
    a.x+=b.x;
    a.y+=b.y;
    a.z+=b.z;
}
void vecAdd(vec &a,vec b,double alpha)
{
    a.x+=b.x*alpha;
    a.y+=b.y*alpha;
    a.z+=b.z*alpha;
}
vec vecAddRet(vec a,vec b,double alpha)
{
    vecAdd(a,b,alpha);
    return a;
}
vec vecSub(vec a,vec b)
{
    vec c;
    c.x=a.x-b.x;
    c.y=a.y-b.y;
    c.z=a.z-b.z;
    return c;
}
vec vecScal(vec a,double s)
{
    a.x*=s;
    a.y*=s;
    a.z*=s;
    return a;
}
vec vecSet(double x,double y,double z)
{
    vec a;
    a.x=x;
    a.y=y;
    a.z=z;
    return a;
}
vec vecSet(double s)
{
    return vecSet(s,s,s);
}

void vecSet(vec &a,double x,double y,double z)
{
    a.x=x;
    a.y=y;
    a.z=z;
}

void vecSet(vec &a,double s)
{
    vecSet(a,s,s,s);
}
void vecRandNeg(vec &a)
{
    double tmp=rand();
    if(tmp>=RAND_MAX/2)
        a.x*=(-1);
    tmp=rand();
    if(tmp>=RAND_MAX/2)
        a.y*=(-1);
    tmp=rand();
    if(tmp>=RAND_MAX/2)
        a.z*=(-1);
}
bool inside(vec a)
{
    if(a.x<0 || a.x>lx)return 0;
    if(a.y<0 || a.y>ly)return 0;
    if(a.z<0 || a.z>lz)return 0;
    return 1;
}
double wrap(double a,double b)
{
    //double ret;
    //if(a<0)
    //    ret=a+double((long long int)(-a/b)+1)*b;
    //else if(a>b)
    //    ret=(long long int)(a)%(long long int)(b-1)-double((long long int)(a))+a;
    //else ret=a;
    //return ret;
    if(a<0)return a+b;
    else if(a>=b)return a-b;
    else return a;
}
void vecWrap(vec &a)
{
    a.x=wrap(a.x,lx);
    a.y=wrap(a.y,ly);
    a.z=wrap(a.z,lz);
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
double len_sq(vec a)
{
    return a.x*a.x+a.y*a.y+a.z*a.z;
}
void histogram()
{
    double h=double(l_hist)/double(d);
    double h_max=0;
    double dv=v_max/double(d-1);

    for(int j=0;j<d;j++)
        x_vel[j]=0;
    for(int i=0;i<n;i++)
        x_vel[int(sqrt(len_sq(vel[i]))/dv)]++;
    for(int j=0;j<d;j++)
        h_max=max(int(h_max),x_vel[j]);

    gclr(win2);
    newcolor(win2,primary.c_str());
    for(int j=0;j<d;j++)
        fillrect(win2,h*double(j)+border2,border2,h,x_vel[j]/h_max*l_hist/2);

    newcolor(win2,bgdark.c_str());
    for(int j=0;j<d;j++)
        drawrect(win2,h*double(j)+border2,border2,h,x_vel[j]/h_max*l_hist/2);

    copylayer(win2,1,0);
}
void redraw()
{
    gclr(win);
    newcolor(win,divider.c_str());
    drawrect(win,double(border),double(border),double(Lx),double(Ly));
    newcolor(win,bglight.c_str());
    fillcirc(win,pos[0].x*mult+double(border),pos[0].y*mult+double(border),5.0,5.0);
    newcolor(win,primary.c_str());
    for(int i=0;i<n;i++)
    {
        if(i)
            fillcirc(win,pos[i].x*mult+double(border),pos[i].y*mult+double(border),5.0,5.0);
    }
    copylayer(win,1,0);
}
void init_graphic()
{
    win=gopen(Lx+border*2,Ly+border*2);
    layer(win,0,1);
    gsetbgcolor(win,bgdark.c_str());

    win2=gopen(l_hist+border2*2,l_hist/2+border2*2);
    gsetbgcolor(win2,bglight.c_str());
    layer(win2,0,1);
}

void energy()
{
    double v_sum=0;
    double vv_sum=0;
    v_max=0;
    for(int i=0;i<n;i++)
    {
        vv_sum+=len_sq(vel[i]);
        v_sum+=sqrt(len_sq(vel[i]));
        v_max=max(v_max,sqrt(len_sq(vel[i])));
    }
    temperature=(1.0/(3.0*n))*vv_sum;
    pressure=density*(vv_sum+vir_sum)/(double(n)*lx*ly*lz);
    kinetic=vv_sum/(2.0*double(n));
    potential/=n;
}

void parameters()
{
    scanf("%d",&n);
    scanf("%d",&N);
    scanf("%lf",&dt);
    scanf("%lf",&density);
    scanf("%lf",&temperature);
    scanf("%d",&Lx);
    scanf("%d",&border);
    scanf("%d",&l_hist);
    scanf("%d",&border2);
    scanf("%d",&d);
    scanf("%d",&h_frames);
    scanf("%d",&wait_time);

    pos.resize(n);
    vel.resize(n);
    acc.resize(n);

    x_vel.resize(d);

    Ly=Lx;
    Lz=Lx;
    mult=pow(double(Lx*Ly*Lz)*density/double(n),1.0/3.0);
    lx=double(Lx)/mult;
    ly=double(Ly)/mult;
    lz=double(Lz)/mult;
}
void colors()
{
    bgdark="#212121";
    bglight="#B2DFDB";
    accent="#00BCD4";
    light="#BDBDBD";
    primary="#009688";
    divider="#BDBDBD";
    darkprimary="#212121";
}

void init()
{
    colors();
    parameters();
    m=int(pow(double(n),1.0/3.0));
    if(m*m*m < n)m++;
    //printf("m=%d n=%d %f\n",m,n,pow(double(n),1.0/3.0));
    
    double dx=lx/double(m);
    double dy=ly/double(m);
    double dz=lz/double(m);

    int u=0;
    int i=0;
    double pi=4*atan(1.0);
    double v_range=sqrt(3*temperature);

    for(int i=0;i<m;i++)
    {
        if(u>=n)break;
        for(int j=0;j<m;j++)
        {
            if(u>=n)break;
            for(int k=0;k<m;k++)
            {
                if(u>=n)break;
                //Random initial positions on a lattice
                pos[u].x=randomRange((double(i)+0.4)*dx,(double(i)+0.6)*dx);
                pos[u].y=randomRange((double(j)+0.4)*dy,(double(j)+0.6)*dy);
                pos[u].z=randomRange((double(k)+0.4)*dz,(double(k)+0.6)*dz);

                // Random initial velocities on a cube
                //vel[u].x=randomRange(0.00001,0.001);
                //vel[u].y=randomRange(0.00001,0.001);
                //vel[u].z=randomRange(0.00001,0.001);
                //vecRandNeg(vel[u]);
                
                //Random initial velocities on a sphere
                double r=randomRange(0.0,v_range);
                //r=v_range;
                double phi=randomRange(0.0,2*pi);
                double tetha=randomRange(0.0,2*pi);
                vel[u].x=r*cos(phi)*sin(tetha);
                vel[u].y=r*sin(phi)*sin(tetha);
                vel[u].z=r*cos(tetha);

                acc[u]=vecSet(0);
                u++;
            }
        }
    }
    energy();
    kinetic0=kinetic;
    potential0=potential;

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
//First Step Leapfrog: v(t+dt/2)=v(t)+(dt/2)a(t)
//                       r(t+dt)=r(t)+dt*v(t+dt/2)
    potential=0;
    double vir_sum=0;
    for(int i=0;i<n;i++)
    {
        vecAdd(vel[i],acc[i],0.5*dt);
        vecAdd(pos[i],vel[i],dt);
        vecWrap(pos[i]);
    }
    for(int i=0;i<n;i++)
        acc[i]=vecSet(0);
    for(int i=0;i<n;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            vec dr_min;
            double dist_min=double(INT_MAX);
            for(int k=0;k<27;k++)
            {
                vec r2=vecAddRet(pos[j],mv[k],1.0);
                vec dr=vecAddRet(pos[i],r2,-1.0);
                double d_ij=len_sq(dr);
                if(dist_min>d_ij)
                {
                    dist_min=d_ij;
                    dr_min=dr;
                }
            }
            double rri=1.0/dist_min;              // (r_ij)^(-2)
            double rri3=rri*rri*rri;              // (r_ij)^(-6)
            double f_ij=48.0*rri3*(rri3-0.5)*rri; // 48*(r_ij^(-13)-0.5*r_ij^(-7)
            vecAdd(acc[i],dr_min,f_ij);
            vecAdd(acc[j],dr_min,-f_ij);

            //Add Potential Energy
            potential+=4.0*rri3*(rri3-1.0)+1.0;
            vir_sum+=f_ij*dist_min;
        }
    }
//Second Step Leapfrog: v(t+dt)=v(t+dt/2)+(dt/2)a(t+dt)
    for(int i=0;i<n;i++)
        vecAdd(vel[i],acc[i],0.5*dt);

//Delete this!! Only for testing
//    for(int i=0;i<n;i++)
//    {
//        vel[i].x=randomRange(25.0,0.0);
//        vel[i].y=randomRange(25.0,0.0);
//        vel[i].z=randomRange(25.0,0.0);
//        vecRandNeg(vel[i]);
//    }
}
void testPos()
{
    for(int i=0;i<n;i++)
        if(!inside(pos[i]))
        {
            printf("outside: ");
            printVec(pos[i]);
        }
}
void ovito()
{
    printf("%d\n",n);
    printf("%d\n",counter);
    for(int i=0;i<n;i++)
        printf("%f %f %f\n",pos[i].x,pos[i].y,pos[i].z);
}
float maxwell()
{
    double sum_vvx=0;
    for(int i=0;i<n;i++)
        sum_vvx+=vel[i].x*vel[i].x;
    return exp(sum_vvx/(n*temperature));
}
void energy_data()
{
    printf("%d %.4f %.4f %.4f %.4f %.4f\n",counter,kinetic,potential,potential+kinetic,temperature,maxwell());
}

int main()
{
    srand(time(NULL));

    init();
    init_graphic();

    counter=0;
    for(int i=0;i<N;i++)
    {
        //printf("temperature: %f\n",temperature);
        potential=0;
        forces();
        energy();
        redraw();
        if(i%h_frames == 0)
            histogram();
        //ovito();
        energy_data();
        counter++;
        if(wait_time)
            msleep(wait_time);
    }
    int key;
    key=ggetch();
    gclose(win);
    return 0;
}
