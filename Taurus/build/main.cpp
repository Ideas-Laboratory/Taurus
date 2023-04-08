#include "layout.h"
#include "Multilevel.h"
#include <sys/time.h>
#include <string>
using namespace std;

int main(int argc, char ** argv){
    //1.定义一个图
    graph g;
    //2.初始化图的信息和结点的初始坐标
    //string s="price_1000.mtx";
    cout<<"11111"<<endl;
    string filename = "bus1138.mtx";
     cout<<"3333"<<endl;
    if (argc > 1) {
		//strcpy(filename, argv[1]);
        filename=argv[1];
	}
    cout<<"222"<<endl;
    string pos="../data/"+filename;
    cout<<"pos="<<pos<<endl;

    g.initgraph(pos,1);
    cout<<"11"<<endl;
    g.initRandomPosition();// pmds  input position
    g.solveDij();
    g.initPivotMDSPosition(200);
    // for(int i=0;i<g.n;i++){
    //     cout<<g.nodes[i].x<<","<<g.nodes[i].y<<endl;
    // }
    //3.自定义力的模型
    vector<force> f;
    f.push_back(force(0,1,1,2));
    f.push_back(force(0,-1,0,1));
    //f.push_back(force(1,40,2,0));
    //f.push_back(force(0,-1,-1,0));
    //4.使用自定义的力求解图布局
    vector<t_force> tf;
    tf.push_back(t_force(1,0.1,1,0,0));
    tf.push_back(t_force(1,0.8,1,0,1));
    tf.push_back(t_force(0,-1,1,0,2));
    // tf.push_back(t_force(1,1,1,0,0));
    // tf.push_back(t_force(1,8,1,0,1));
    // tf.push_back(t_force(0,-10,1,0,2));
    struct timeval start,end;
    gettimeofday(&start,NULL);
    //FMMMLayout(g);
   // FDPLayout(g);
    Layout(g,f);
    //t_Layout(g,tf);
    gettimeofday(&end,NULL);
    float time_use=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);//微秒
    printf("time_use is %.3fs\n",time_use/1000000);
   //FMMMLayout fmmm;
    //fmmm.call(g);

    //f.push_back(force(1,1,1,1));
    //f.push_back(force(3,1,1,-1));
    //f.push_back(force(0,-1,-1,0));
    //SMLayout(g);
    

    double edgelen = 0;
    int edgenum = 0;
    float sum_ddij = 0, sum_d = 0, scale = 1;
    for (int i = 0; i < g.nodes.size(); i++) {
        for (int j = 0; j < g.nodes.size(); j++) {
            if(g.sp[i][j]==0) continue;
            float disx = g.nodes[i].x - g.nodes[j].x, disy = g.nodes[i].y - g.nodes[j].y;
            sum_ddij += (sqrt(disx * disx + disy * disy) * g.sp[i][j])/(g.sp[i][j]*g.sp[i][j]);
            sum_d += (disx * disx + disy * disy)/(g.sp[i][j]*g.sp[i][j]);
        }
    }
    scale = sum_ddij / sum_d;
    for (int i = 0; i < g.nodes.size(); i++) {
        for (int j = 0; j < i; j++) {
            if(g.sp[i][j]==0) continue;
            double d = sqrt(pow(scale*(g.nodes[i].x - g.nodes[j].x), 2) + pow(scale*(g.nodes[i].y - g.nodes[j].y), 2));
            double dis = pow(d - g.sp[i][j], 2);
            edgelen += dis / pow(g.sp[i][j], 2);
            edgenum++;
        }
    }
   // string output= "fig/" +filename;
    g.drawSVG(filename,"t-FDP",1);
   // g.drawSVG(output,"t-FDP",1);
    cout << "stress:" <<setiosflags(ios::fixed)<<setprecision(3)<< (edgelen / edgenum) << endl;
    cout<<"NP1="<<testNP3(g.nodes, g.links, g.sp,1)<<endl;
    cout<<"NP2="<<testNP3(g.nodes, g.links, g.sp,2)<<endl;

}