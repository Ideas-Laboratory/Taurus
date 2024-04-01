#include "metric.h"
double testSE(std::vector<Node> &nodes, std::vector<Link> &links, std::vector<std::vector<double>> &sp)
{
    double edgelen = 0;
    int edgenum = 0;
    float sum_ddij = 0, sum_d = 0, scale = 1;
    for (int i = 0; i < nodes.size(); i++)
    {
        for (int j = 0; j < nodes.size(); j++)
        {
            if (sp[i][j] == 0)
                continue;
            float disx = nodes[i].x - nodes[j].x, disy = nodes[i].y - nodes[j].y;
            sum_ddij += (sqrt(disx * disx + disy * disy) * sp[i][j]) / (sp[i][j] * sp[i][j]);
            sum_d += (disx * disx + disy * disy) / (sp[i][j] * sp[i][j]);
        }
    }
    scale = sum_ddij / sum_d;
    for (int i = 0; i < nodes.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (sp[i][j] == 0)
                continue;
            double d = sqrt(pow(scale * (nodes[i].x - nodes[j].x), 2) + pow(scale * (nodes[i].y - nodes[j].y), 2));
            double dis = pow(d - sp[i][j], 2);
            edgelen += dis / pow(sp[i][j], 2);
            edgenum++;
        }
    }
    return edgelen / edgenum;
}
float testNP(vector<Node> &nodes,vector<Link> &links,vector<vector<double>> &sp,int neigh){
    double NP = 0;
    int leafnum = 0;
    int number=nodes.size();
    int countt=0;
    for (int i = 0; i < number; i++) {
        vector<int> dijknn;
        vector<double> dknn;
        for (int j = 0; j < number; j++) {
            if (j != i) {
                if (sp[i][j] <= neigh&&sp[i][j]>0) {
                    dijknn.push_back(j);
                }
            }
        }
        for (int j = 0; j < dijknn.size(); j++) {
            dknn.push_back(1000000000000);
        }
        for (int j = 0; j < number; j++) {
            if (j != i) {
                double curdis = sqrt(pow(nodes[i].x - nodes[j].x, 2) + pow(nodes[i].y - nodes[j].y, 2));
                for (int k = 0; k < dijknn.size(); k++) {
                    if (curdis < dknn[k]) {
                        for (int l = dijknn.size() - 1; l > k; l--) {
                            dknn[l] = dknn[l - 1];
                        }
                        dknn[k] = curdis;
                        break;
                    }
                }
            }
        }
        int intersect = 0;
        int uni = 0;
        for (int j = 0; j < dijknn.size(); j++) {
            double curdis = sqrt(pow(nodes[i].x - nodes[dijknn[j]].x, 2) + pow(nodes[i].y - nodes[dijknn[j]].y, 2));
            if (curdis <= dknn[dknn.size() - 1]) {
                intersect++;
            }
            else {
                uni += 2;
            }
        }
        if (dijknn.size() == 0) { NP += 0; leafnum += 1; }
        else { NP += (double)intersect / (uni + intersect); }
        countt+=intersect;
    }
    cout << "NP" <<neigh<<"="<<setiosflags(ios::fixed)<<setprecision(3)<<NP / (double)number << endl;
    return NP/(double)number;
}



double testCD(vector<Node> &nodes,vector<Link> &links,vector<int> &flag,int clusternum){
    //cout<<"clusterdata"<<nodes.size()<<" "<<links.size()<<" "<<flag.size()<<endl;
    float maxx=0,minx=100000000000,maxy=0,miny=10000000000;

    for(int i=0;i<nodes.size();i++){
        if(nodes[i].x>maxx) maxx = nodes[i].x;
        if(nodes[i].x<minx) minx = nodes[i].x;
        if(nodes[i].y>maxy) maxy = nodes[i].y;
        if(nodes[i].y<miny) miny = nodes[i].y;
    }
    float lenth=max(maxx-minx,maxy-miny);
    for(int i=0;i<nodes.size();i++){
        nodes[i].x-=minx;
        nodes[i].y-=miny;
    }
    for(int i=0;i<nodes.size();i++){
        nodes[i].x/=lenth;
        nodes[i].y/=lenth;
//        nodes[i].x*=10;
//        nodes[i].y*=10;
    }
    vector<vector<Node>> all_clusters;
    for(int i=0;i<clusternum;i++){
        vector<Node> cluster;
        all_clusters.push_back(cluster);
    }
    for(int i=0;i<nodes.size();i++){
        all_clusters[flag[i]].push_back(nodes[i]);
    }
    float average_inter=0;
    int inter_num=0;

    for(int i =0;i<nodes.size();i++){
        float mindis=1;
        for(int j =0;j<nodes.size();j++){
            if(flag[i]==flag[j])continue;
            float dis= sqrt((nodes[i].x-nodes[j].x)*(nodes[i].x-nodes[j].x)+(nodes[i].y-nodes[j].y)*(nodes[i].y-nodes[j].y));
            if(dis<mindis){
                mindis=dis;
            }

        }
        average_inter+=mindis;
        inter_num++;

    }
    average_inter/=inter_num;
    cout << "CD" <<"="<<setiosflags(ios::fixed)<<setprecision(3)<<average_inter<< endl;
    return average_inter;
}
double testCE(vector<Node> &nodes,vector<Link> &links,vector<int> &flag,int clusternum){
    //cout<<"clusterdata"<<nodes.size()<<" "<<links.size()<<" "<<flag.size()<<endl;
    float maxx=0,minx=100000000000,maxy=0,miny=10000000000;

    for(int i=0;i<nodes.size();i++){
        if(nodes[i].x>maxx) maxx = nodes[i].x;
        if(nodes[i].x<minx) minx = nodes[i].x;
        if(nodes[i].y>maxy) maxy = nodes[i].y;
        if(nodes[i].y<miny) miny = nodes[i].y;
    }
    float lenth=max(maxx-minx,maxy-miny);
    for(int i=0;i<nodes.size();i++){
        nodes[i].x-=minx;
        nodes[i].y-=miny;
    }
    for(int i=0;i<nodes.size();i++){
        nodes[i].x/=lenth;
        nodes[i].y/=lenth;
//        nodes[i].x*=10;
//        nodes[i].y*=10;
    }
    vector<vector<Node>> all_clusters;
    for(int i=0;i<clusternum;i++){
        vector<Node> cluster;
        all_clusters.push_back(cluster);
    }
    for(int i=0;i<nodes.size();i++){
        all_clusters[flag[i]].push_back(nodes[i]);
    }
    float average_inter=0;
    for(int k=0;k<clusternum;k++){
        float clusterdis=0;
        int clusternum=0;
        for(int i=0;i<all_clusters[k].size();i++){
            for(int j=0;j<i;j++){
                float dis=sqrt((all_clusters[k][i].x-all_clusters[k][j].x)*(all_clusters[k][i].x-all_clusters[k][j].x)+(all_clusters[k][i].y-all_clusters[k][j].y)*(all_clusters[k][i].y-all_clusters[k][j].y));
                clusterdis+=dis;
                clusternum++;
            }
        }
        clusterdis/=clusternum;
        average_inter+=clusterdis;
    }
    average_inter/=clusternum;
    cout << "CE" <<"="<<setiosflags(ios::fixed)<<setprecision(3)<<average_inter<< endl;
    return average_inter;
}


int judgeCross(Node x,Node y,Node z,Node h){
    float m=(x.x-z.x)*(h.y-z.y)-(x.y-z.y)*(h.x-z.x);
    float n=(y.x-z.x)*(h.y-z.y)-(y.y-z.y)*(h.x-z.x);
    if(m*n<=0) return 1;
    return 0;
}
float multi(Node p1,Node p2,Node p0){
    return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
Node intersect(Node a,Node b,Node c,Node d){
    Node p;
    p.x=(multi(a,d,c)*b.x-multi(b,d,c)*a.x)/(multi(a,d,c)-multi(b,d,c));
    p.y=(multi(a,d,c)*b.y-multi(b,d,c)*a.y)/(multi(a,d,c)-multi(b,d,c));
    return p;
}
double testCL(vector<Node> &nodes,vector<Link> &links){
    int V=nodes.size(),E=links.size();
    float n_all=E*(E-1)/2;
    float n_imp=0,c_max=0,cn=0;
    float crosslessness=1;
    //compute degree
    vector<int> count(V,0);
    for(int i=0;i<E;i++){
        count[links[i].source]++;
        count[links[i].target]++;
    }

    for(int i=0;i<V;i++){
        n_imp+=count[i]*(count[i]-1);
    }
    n_imp/=2;
    c_max=n_all-n_imp;
    if(c_max==0)cout<<"cross error"<<n_all<<" "<<E<<endl;

//#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<E;i++){
        int src=links[i].source,tgt=links[i].target;
        if(src==tgt) continue;

        float src_x=nodes[src].x,src_y=nodes[src].y;
        float tgt_x=nodes[tgt].x,tgt_y=nodes[tgt].y;
        for(int j=i+1;j<E;j++){
            int temp_src=links[j].source,temp_tgt=links[j].target;
            if(src==temp_src||tgt==temp_src||src==temp_tgt||tgt==temp_tgt) continue;
            if(temp_tgt==temp_src) continue;
            //compute if cross
            Node a,b,c,d;
            a.x=src_x,a.y=src_y;
            b.x=tgt_x,b.y=tgt_y;
            c.x=nodes[temp_src].x,c.y=nodes[temp_src].y;
            d.x=nodes[temp_tgt].x,d.y=nodes[temp_tgt].y;

            int ans1=judgeCross(a,b,c,d);
            int ans2=judgeCross(c,d,a,b);

            if(ans1==1&&ans2==1){
                cn++;
            }

        }
    }
    crosslessness=1-sqrtf(cn/c_max);
    cout<<"cross="<<crosslessness<<endl;
    return crosslessness;
}
double testMA(vector<Node> &nodes,vector<Link> &links){
    int V=nodes.size(),E=links.size();
    double min_angle=0,sum_angle=0;

    vector<Edge> e[V];
    int st, end, weight;

    for (int i=0; i < E; i++)
    {
        //cin>>st>>end>>weight;
        st=links[i].source;
        end=links[i].target;
        Edge tmp;
        //tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);    //vector
        st=links[i].target;
        end=links[i].source;
        //tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    for(int i=0;i<V;i++){
        if(e[i].size()<2) continue;
        double theta=(double)360/(double)e[i].size();
        double min_v_angle=360;
        //cout<<count[i]<<" "<<e[i].size()<<" ";
        for(int j=0;j<e[i].size()-1;j++){
            int v0=i;
            int v1=e[i][j].end;
            double v1_x=nodes[v0].x-nodes[v1].x;
            double v1_y=nodes[v0].y-nodes[v1].y;
            for(int k=j+1;k<e[i].size();k++){
                int v2=e[i][k].end;
                //compute angel
                double v2_x=nodes[v0].x-nodes[v2].x;
                double v2_y=nodes[v0].y-nodes[v2].y;
                if((sqrt(v1_x*v1_x+v1_y*v1_y)*sqrt(v2_x*v2_x+v2_y*v2_y))==0) continue;
                double angle=((v1_x*v2_x)+(v1_y*v2_y))/(sqrt(v1_x*v1_x+v1_y*v1_y)*sqrt(v2_x*v2_x+v2_y*v2_y));
                if(angle<=1&&angle>=-1){
                    angle=acos(angle)*180/PI;
                }

                if(min_v_angle>angle) min_v_angle=angle;
            }
        }

        double temp_angle=abs((theta-min_v_angle)/theta);
        sum_angle+=temp_angle;
        //cout<<"theta="<<theta<<"min angle at v="<<abs((theta-min_v_angle)/theta)<<"sum_angle="<<sum_angle<<"        ";
    }
    //cout<<endl;
    //cout<<"sum_angle="<<sum_angle<<endl;
    min_angle=1-sum_angle/(double)V;
    cout<<"min_angle="<< min_angle<<endl;
    return min_angle;
}
