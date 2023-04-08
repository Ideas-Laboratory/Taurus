#include "solver.h"

int str2int(std::string s) {
    std::stringstream stream(s);
    int i;
    stream >> i;
    stream.clear();
    return i;
}
float str2float(std::string s) {
    std::stringstream stream(s);
    float i;
    stream >> i;
    stream.clear();
    return i;
}


std::vector<int> split(const std::string& str, const std::string& delim) {
    std::vector<int> res;
    if ("" == str) return res;
    char * strs = new char[str.length() + 1];
    strcpy(strs, str.c_str());

    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while (p) {
        std::string s = p;
        int y = atoi(s.c_str());
        res.push_back(y);
        p = strtok(NULL, d);
    }
    return res;
}
std::vector<std::string> split_string(const std::string& str, const std::string& delim) {
    std::vector<std::string> res;
    if ("" == str) return res;
    char * strs = new char[str.length() + 1];
    strcpy(strs, str.c_str());

    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while (p) {
        std::string s = p;
       // int y = atoi(s.c_str());
        res.push_back(s);
        p = strtok(NULL, d);
    }
    return res;
}


graph::graph() {
    n=0;m=0;
}

/*void graph::readGraph(std::string filename){
    ifstream input(filename);
    std::string s;
    std::string b="break";
    while (getline(input, s)) {
        if (s.at(0) != '%') {
            break;
        }
    }
    std::cout<<s<<std::endl;
    std::vector<int> di = split(s, " ");
    if(di.size()==3){
        n = di[1];
        m = di[2];

        //std::vector<Link> links = std::vector<Link>(m);
        for(int i=0;i<m;i++){
            if(!getline(input,s)) std::cout<<"file read error!"<<std::endl;
            std::vector<int> di = split(s, " ");
            Link temp;
            temp.source=di[0]-1;
            temp.target=di[1]-1;
            links.push_back(temp);
        }
    }else{
        int max=0;
        int eps=0;
        *//*if(di[0]==0||di[1]==0){
            eps=0;
        }*//*
        Link temp;
        temp.source=di[0]-eps;
        temp.target=di[1]-eps;
        if(max<di[0]) max=di[0]+1-eps;
        if(max<di[1]) max=di[1]+1-eps;
        links.push_back(temp);
        while (getline(input, s)) {
            if(s==b){break;}
            std::vector<int> di1 = split(s, " ");
            Link temp;
            temp.source=di1[0]-eps;
            temp.target=di1[1]-eps;
            links.push_back(temp);
            if(max<di1[0]) max=di1[0]+1-eps;
            if(max<di1[1]) max=di1[1]+1-eps;
        }
        m=links.size();
        n=max;
    }
    std::cout<<n<<" "<<m<<std::endl;
}*/

void graph::initgraph(std::string filename) {
    std::ifstream input(filename);
    std::string s;
    std::string b="break";
    while (getline(input, s)) {
        if (s.at(0) != '%') {
            break;
        }
    }
    std::vector<int> I,J;

    std::vector<int> di = split(s, " ");
    if(di.size()==3){
        n = di[1];
        m = di[2];
    }else{
        I.push_back(di[0]);
        J.push_back(di[1]);
    }
    while(getline(input,s)){
        di=split(s, " ");
        I.push_back(di[0]);
        J.push_back(di[1]);
    }
    int max_i=*max_element(I.begin(),I.end());
    int max_j=*max_element(J.begin(),J.end());
    int min_i=*min_element(I.begin(),I.end());
    int min_j=*min_element(J.begin(),J.end());
    int max_ind=max_i>max_j?max_i:max_j;
    int min_ind=min_i<min_j?min_i:min_j;
    n=max_ind-min_ind+1;
    m=I.size();
    std::vector<std::set<int> > undirected(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if edge not seen
        {
            undirected[i].insert(j);
            undirected[j].insert(i);
            Link temp;
            temp.source=i-min_ind;
            temp.target=j-min_ind;
            links.push_back(temp);
        }
    }

    m=links.size();
    nodes.resize(n);
    std::cout<<"build undirected graph n="<<n<<" m="<<m<<std::endl;
}
void graph::initgraph(std::string filename,bool stringid) {
    std::ifstream input(filename);
    std::string s;
    std::string b="break";
    while (getline(input, s)) {
        if (s.at(0) != '%'&&s.at(0) != '#') {
            break;
        }
    }
   // getline(input, s);
 
    //std::vector<int> quasi_node(46);
    //std::vector<int> fflag={2741, 12781, 3372, 6610, 4164, 4511, 2212, 6830, 19961, 15003, 15659, 20635, 25758, 14807, 570, 24955, 4513, 6179, 14540, 17655, 20108, 19423, 25346, 22691, 12496, 22887, 7956, 1653, 2952, 45, 12851, 18894, 23293, 9785, 8879, 46, 20562, 21012, 12365, 21508, 17692, 11472, 773, 21847, 11241, 21281};
   // std::vector<int> flag_id(fflag.size(),0);
   // std::vector<int> fflag={0,1,10,12,2,3,15,18,19,21};
    std::vector<int> I,J;
    std::vector<int> di = split(s, " ");
    std::map<int,int> indexById;
    int id=0;
    if(di.size()==3){
        n = di[1];
        m = di[2];
    }else{
        I.push_back(di[0]);
        J.push_back(di[1]);
        if(indexById.count(di[0])==0){
            indexById.insert(std::pair<int,int>(di[0],id));
            id++;
        }
        if(indexById.count(di[1])==0){
            indexById.insert(std::pair<int,int>(di[1],id));
            id++;
        }
        std::cout<<di[0]<<","<<di[1]<<std::endl;
    }
    while(getline(input,s)){
        di=split(s, " ");
        I.push_back(di[0]);
        J.push_back(di[1]);
        if(indexById.count(di[0])==0){
            indexById.insert(std::pair<int,int>(di[0],id));
            id++;
        }
        if(indexById.count(di[1])==0){
            indexById.insert(std::pair<int,int>(di[1],id));
            id++;
        }
    }
    std::cout<<"id="<<id<<std::endl;
    std::vector<int> id_I,id_J;
    for(int i=0;i<I.size();i++){
        int _i=indexById[I[i]];
        int _j=indexById[J[i]];
        id_I.push_back(_i);
        id_J.push_back(_j);
    }
    int max_i=*max_element(id_I.begin(),id_I.end());
    int max_j=*max_element(id_J.begin(),id_J.end());
    int min_i=*min_element(id_I.begin(),id_I.end());
    int min_j=*min_element(id_J.begin(),id_J.end());
    int max_ind=max_i>max_j?max_i:max_j;
    int min_ind=min_i<min_j?min_i:min_j;
    n=max_ind-min_ind+1;
    m=id_I.size();
    std::vector<std::set<int> > undirected(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = id_I[ij], j = id_J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if edge not seen
        {
            undirected[i].insert(j);
            undirected[j].insert(i);
            Link temp;
            temp.source=i-min_ind;
            temp.target=j-min_ind;
            links.push_back(temp);
        }
    }

    m=links.size();
    nodes.resize(n);
    // std::cout<<"{";
    // for(int i=0;i<fflag.size();i++){
    //     std::cout<<indexById[fflag[i]]<<",";
    // }
    // std::cout<<"}";
    std::cout<<"build undirected graph n="<<n<<" m="<<m<<std::endl;
    max_connect();
}
void graph::max_connect(){
    std::vector<int> ind(n,-1);
    std::vector<int> count_com(n,0);
    std::vector<Edge> e[n];
    int st, end, weight;
   // std::vector<int> fflag={265,282,267,72,269,270,264,274,291,159,286,294,303,285,261,301,271,273,284,288,292,100,302,103,280,299,275,263,266,259,283,290,300,77,276,260,293,101,279,296,289,278,262,297,277,295};
  // std::vector<int> fflag={0,1,10,12,2,3,15,18,19,21};
    // std::vector<int> flag_ind(n,-1);
    // for(int i=0;i<fflag.size();i++){
    //     flag_ind[fflag[i]]=1;
    // }
    for (int i=0; i < m; i++) {
        st = links[i].source;
        end = links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
        st = links[i].target;
        end = links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    int count=0;
    for(int i=0;i<n;i++){
        if(ind[i]==-1){
            ind[i]=count;
            //cout<<"indi=="<<ind[i]<<endl;
            BFS(n,m,e,i,ind,count_com);
            count++;
        }
    }
 if(count>1){
    std::vector<int> flag_set;
 //if(0){
        std::cout<<"The graph is disconnected! I am taking the largest component!"<<std::endl;

    int max_com=0,max=0;
    //int count_con=0;
    for(int i=0;i<count;i++){
        int temp=count_com[i];
        if(temp>max){max=temp;max_com=i;}
    }
    std::cout<<"The component has="<<count_com[max_com]<<" vertex!"<<std::endl;
    max_com=ind[0];
    count=0;
    std::map<int,int> map_ind;
    for(int i=0;i<n;i++){
        if(ind[i]!=max_com){
            count++;
            map_ind.insert(std::pair<int,int>(i,-1));
        }else{
            map_ind.insert(std::pair<int,int>(i,i-count));
            //if(flag_ind[i]==1){flag_set.push_back(i-count);}
        }
    }
    //cout<<"count=="<<count<<endl;
    n-=count;
    int m_count=0;
    for(int i=0;i<links.size();i++){
        int src=links[i].source,tgt=links[i].target;
        if(map_ind[src]==-1||map_ind[tgt]==-1){m_count++;}
        if(map_ind[src]!=-1){
            links[i-m_count].source=map_ind[src];
        }
        if(map_ind[tgt]!=-1){
            links[i-m_count].target=map_ind[tgt];
        }

    }
   // std::cout<<"{";
    // for(int i=0;i<flag_set.size();i++){
    //     std::cout<<flag_set[i]<<",";
    // }
   
    // int fcount=0;
    // for(int i=0;i<flag.size();i++){
    //     if(map_ind[i]==-1){fcount++;}else{
    //         flag[i-fcount]=flag[i];
    //     }
    // }

    m-=m_count;
    
    std::string filename1="fGSE.mtx";
    freopen(filename1.c_str(),"w",stdout);
    for(int i=0;i<m;i++){
        std::cout<<links[i].source<<" "<<links[i].target<<std::endl;
    }
  //  std::cout<<"break"<<std::endl;
    // for(int i=0;i<n;i++){
    //     cout<<flag[i]<<endl;
    // }
   // cin.clear();
   // cout.clear();
    freopen("/dev/tty", "w", stdout);
    //  freopen("CON", "r", stdin);
    //  freopen("CON", "w", stdout);

}
}

void graph::readGraphMtx(std::string filename){
    std::ifstream input(filename);
    std::string s;
    std::string b="break";
    int src,tgt,num_n=0,num_m=0;
    while(getline(input, s)){
        if(s==b){break;}
        num_m++;
        std::vector<int> di = split(s, " ");
        Link temp;
        src=di[0];
        tgt=di[1];
        temp.source=src;
        temp.target=tgt;
        links.push_back(temp);
        num_n=std::max(num_n,src);
        num_n=std::max(num_n,tgt);
    }

    n=num_n+1;
    m=num_m;

    std::vector<int> ind(n,-1);
    std::vector<int> count_com(n,0);

    //max component
    //std::vector<Edge> e[n];
    std::vector<Edge> *e=new std::vector<Edge>[n];
    int st, end, weight;
    for (int i=0; i < m; i++) {
        st = links[i].source;
        end = links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
        st = links[i].target;
        end = links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    int count=0;
    for(int i=0;i<n;i++){
        if(ind[i]==-1){
            ind[i]=count;
            //std::cout<<"indi=="<<ind[i]<<std::endl;
            BFS(n,m,e,i,ind,count_com);
            count++;
        }
    }
    if(count>1){
        //std::cout<<"The graph is disconnected! I am taking the largest component!"<<std::endl;

        int max_com=0,max=0;
        for(int i=0;i<count;i++){
            int temp=count_com[i];
            if(temp>max){max=temp;max_com=i;}
        }
        //std::cout<<"The component has="<<count_com[max_com]<<" vertex!"<<std::endl;

        count=0;
        std::map<int,int> map_ind;
        for(int i=0;i<n;i++){
            if(ind[i]!=max_com){
                count++;
                map_ind.insert(std::pair<int,int>(i,-1));
            }else{
                map_ind.insert(std::pair<int,int>(i,i-count));
            }
        }
        //std::cout<<"count=="<<count<<std::endl;
        n-=count;
        int m_count=0;
        for(int i=0;i<links.size();i++){
            int src=links[i].source,tgt=links[i].target;
            if(map_ind[src]==-1||map_ind[tgt]==-1){m_count++;}
            if(map_ind[src]!=-1){
                links[i-m_count].source=map_ind[src];
            }
            if(map_ind[tgt]!=-1){
                links[i-m_count].target=map_ind[tgt];
            }
        }

        m-=m_count;

    }

}

void graph::initPosition() {
    //nodes[i].x nodes[i].y
    if(nodes.size()!=n){nodes.resize(n);}
    for (int i = 0; i <this->n; i++) {
        double radius = 10.0 * sqrt(i), angle = i * 3.141592653589793f * (3.0 - sqrt(5.0));
        Node temp;
        temp.id=i;
        temp.x = radius * cos(angle);
        temp.y = radius * sin(angle);
        nodes[i]=temp;
    }
    for(int i=0;i<this->m;i++){
        nodes[links[i].source].deg++;
        nodes[links[i].target].deg++;
    }
}
void graph::initRandomPosition() {
    if(nodes.size()!=n){nodes.resize(n);}
    for (int i = 0; i <this->n; i++) {
        double radius = 10.0 * sqrt(i), angle = i * 3.141592653589793f * (3.0 - sqrt(5.0));
        Node temp;
        temp.id=i;
        temp.x = radius * cos(angle);
        temp.y = radius * sin(angle);
        nodes[i]=temp;
    }
    for(int i=0;i<this->m;i++){
        nodes[links[i].source].deg++;
        nodes[links[i].target].deg++;
    }
}
void graph::initPivotMDSPosition(int p){

    //if(nodes.size()!=n){nodes.resize(n);}
    if(p>n) p=n;
    int I[m];
    int J[m];
    for (int i = 0; i < m; i++) {
        I[i]=links[i].source;
        J[i]=links[i].target;
    }
    std::vector<std::vector<int> > g1 = build_graph_unweighted(n, m, I, J);
    std::vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 48);
    std::set<int> pset;
    for(int i=0;i<closest_pivots.size();i++){
        pset.insert(closest_pivots[i]);
        //std::cout<<closest_pivots[i]<<",";
    }
    std::vector<int> pivot;
    std::cout<<std::endl<<"size="<<pset.size()<<std::endl;
    for(std::set<int>::iterator it=pset.begin();it!=pset.end();it++){
        pivot.push_back(*it);
    }

    std::cout<<"compute pivotMDS"<<std::endl;
    PivotMDS(nodes,pivot.size(),n,sp,pivot);
    float maxx=0,minx=100000000000,maxy=0,miny=10000000000;

    for(int i=0;i<nodes.size();i++){
        if(nodes[i].x>maxx) maxx = nodes[i].x;
        if(nodes[i].x<minx) minx = nodes[i].x;
        if(nodes[i].y>maxy) maxy = nodes[i].y;
        if(nodes[i].y<miny) miny = nodes[i].y;
    }
    float lenth=max(maxx-minx,maxy-miny);
    float lenth_x=maxx-minx,lenth_y=maxy-miny;
    for(int i=0;i<nodes.size();i++){
        nodes[i].x-=minx;
        nodes[i].y-=miny;
    }
    for(int i=0;i<nodes.size();i++){
        nodes[i].x/=lenth_x;
        nodes[i].y/=lenth_y;
       // nodes[i].x*=100.0f;
        //nodes[i].y*=100.0f;
    }
}
void graph::initNode(int n) {
    nodes.resize(n);
}
void graph::initEdge(int m) {
    links.resize(m);
}


void graph::solveDij() {

    std::vector<Edge> e[n];
    int st, end, weight;
    for (int i=0; i < this->m; i++) {

        st = links[i].source;
        end = links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
        st = links[i].target;
        end = links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
//#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<n;i++){
        std::vector<float> dist(n);
        Dij(n,m,e,i,dist);
        std::vector<double> isp(n);
        for(int j=0;j<n;j++){
            if(dist[j]==MAX){
                isp[j]=0;
            } else{
                isp[j]=dist[j];
            }

        }
        sp.push_back(isp);
    }

}


void graph::append_E_range(para pa) {
    //
    //consSGD.insert();

    for(int i=0;i<m;i++){
        int src=links[i].source,tgt=links[i].target;
        if(src!=tgt){
            float w;
            if(sp[src][tgt]!=0){
                w= 1.0f / (sp[src][tgt] * sp[src][tgt]);
            }else{
                w = 0.0;
            }
            constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));
            if(pivots.size()){
                std::cout<<"pivot"<<std::endl;
                constraints.push_back(constraint(tgt,src,sp[src][tgt],w,pa));//np
            }

        }
    }
}
void graph::append_E_range(force pa) {
    //
    //consSGD.insert();

    for(int i=0;i<m;i++){
        int src=links[i].source,tgt=links[i].target;
        if(src!=tgt){
            float w;
            if(sp[src][tgt]!=0){
                w= 1.0f / (sp[src][tgt] * sp[src][tgt]);
            }else{
                w = 0.0;
            }
            if(src>tgt){
                //constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));
                int temp=src;
                src=tgt;
                tgt=temp;
            }
            if(append_cons[src][tgt].flag){
                append_cons[src][tgt].force_m.push_back(pa);
            }else{
                append_cons[src][tgt].flag=true;
                append_cons[src][tgt].d=sp[src][tgt];
                append_cons[src][tgt].w=w;
                append_cons[src][tgt].force_m.push_back(pa);
            }
           // auto it = std::set_sgd.find(constraint_sgd(src,tgt));
            // if(it!=std::set_sgd.end()){
            //     it->force_sgd.push_back(pa);
            // }else{
            //     std::vector<force> f;
            //     f.push_back(pa);
            //     std::set_sgd.insert(constraint_sgd(src,tgt,sp[src][tgt],w,f));
            // }

            /*if(pivots.size()){
                std::cout<<"pivot"<<std::endl;
                constraints.push_back(constraint(tgt,src,sp[src][tgt],w,pa));//np
            }*/

        }
    }
}
void graph::append_E_range(t_force pa) {
    //
    //consSGD.insert();

    for(int i=0;i<m;i++){
        int src=links[i].source,tgt=links[i].target;
        if(src!=tgt){
            float w;
            if(sp[src][tgt]!=0){
                w= 1.0f / (sp[src][tgt] * sp[src][tgt]);
            }else{
                w = 0.0;
            }
            if(src>tgt){
                int temp=src;
                src=tgt;
                tgt=temp;
            }
            if(append_cons[src][tgt].flag){
                append_cons[src][tgt].t_force_m.push_back(pa);
            }else{
                append_cons[src][tgt].flag=true;
                append_cons[src][tgt].d=sp[src][tgt];
                append_cons[src][tgt].w=w;
                append_cons[src][tgt].t_force_m.push_back(pa);
            }

        }
    }
}

void graph::append_N_range(para pa) {
    //
    //consSGD.insert()
    for(int i=0;i<n;i++){
        for(int j=0;j<i;j++){
            if(i!=j){
                float w;
                if(sp[i][j]!=0){
                    w = 1.0f / (sp[i][j] * sp[i][j]);
                }else{
                    w = 0.0;
                }
                constraints.push_back(constraint(i,j,sp[i][j],w,pa));
            }

        }
    }

}
void graph::append_N_range(force pa) {
    //
    //consSGD.insert()
    for(int i=0;i<n;i++){
        for(int j=0;j<i;j++){
            if(i!=j){
                float w;
                if(sp[i][j]!=0){
                    w = 1.0f / (sp[i][j] * sp[i][j]);
                }else{
                    w = 0.0;
                }
                int tgt=j,src=i;
                 if(src>tgt){
                //constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));
                int temp=src;
                src=tgt;
                tgt=temp;
                }
                if(append_cons[src][tgt].flag){
                     append_cons[src][tgt].force_m.push_back(pa);
                }else{
                    append_cons[src][tgt].flag=true;
                    append_cons[src][tgt].d=sp[src][tgt];
                    append_cons[src][tgt].w=w;
                    append_cons[src][tgt].force_m.push_back(pa);
                }
                //constraints.push_back(constraint(i,j,sp[i][j],w,pa));
                // auto it = std::set_sgd.find(constraint_sgd(i,j));
                // if(it!=std::set_sgd.end()){
                //     it->force_sgd.push_back(pa);
                // }else{
                //     std::vector<force> f;
                //     f.push_back(pa);
                //     std::set_sgd.insert(constraint_sgd(i,j,sp[i][j],w,f));
                // }
            }

        }
    }

}
void graph::append_N_range(t_force pa) {
    //
    //consSGD.insert()
    for(int i=0;i<n;i++){
        for(int j=0;j<i;j++){
            if(i!=j){
                float w;
                if(sp[i][j]!=0){
                    w = 1.0f / (sp[i][j] * sp[i][j]);
                }else{
                    w = 0.0;
                }
                int tgt=j,src=i;
                 if(src>tgt){
                //constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));
                int temp=src;
                src=tgt;
                tgt=temp;
                }
                if(append_cons[src][tgt].flag){
                     append_cons[src][tgt].t_force_m.push_back(pa);
                }else{
                    append_cons[src][tgt].flag=true;
                    append_cons[src][tgt].d=sp[src][tgt];
                    append_cons[src][tgt].w=w;
                    append_cons[src][tgt].t_force_m.push_back(pa);
                }
            }

        }
    }

}
void graph::append_range(std::vector<Link> range,para pa){
    for(int i=0;i<range.size();i++){
        int src=range[i].source,tgt=range[i].target;
        if(src!=tgt){
            float w;
            if(sp[src][tgt]!=0){
                w= 1.0f / (sp[src][tgt] * sp[src][tgt]);
            }else{
                w = 0.0;
            }
            constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));


        }
    }
}
void graph::append_range(std::vector<Link> range,force pa){
    for(int i=0;i<range.size();i++){
        int src=range[i].source,tgt=range[i].target;
        if(src!=tgt){
            float w;
            if(sp[src][tgt]!=0){
                w= 1.0f / (sp[src][tgt] * sp[src][tgt]);
            }else{
                w = 0.0;
            }
            // if(src<tgt){
            //     int temp=src;
            //     src=tgt;
            //     tgt=temp;
            // }
             if(src>tgt){
                //constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));
                int temp=src;
                src=tgt;
                tgt=temp;
            }
            if(append_cons[src][tgt].flag){
                append_cons[src][tgt].force_m.push_back(pa);
            }else{
                append_cons[src][tgt].flag=true;
                append_cons[src][tgt].d=sp[src][tgt];
                append_cons[src][tgt].w=w;
                append_cons[src][tgt].force_m.push_back(pa);
            }
            // auto it = std::set_sgd.find(constraint_sgd(src,tgt));
            // if(it!=std::set_sgd.end()){
            //     it->force_sgd.push_back(pa);
            // }else{
            //     std::vector<force> f;
            //     f.push_back(pa);
            //     std::set_sgd.insert(constraint_sgd(src,tgt,sp[src][tgt],w,f));
            // }
        }
    }
}
void graph::append_Np_range(para pa) {

    //consSGD.insert()
    std::vector<constraint> constraintnp;
    int p=100;
    int I[m];
    int J[m];
    for (int i = 0; i < m; i++) {
        I[i]=links[i].source;
        J[i]=links[i].target;
    }
    pivots.resize(n);
    for(int i=0;i<pivots.size();i++){
        pivots[i]=0;
    }
    std::vector<std::vector<int> > g1 = build_graph_unweighted(n, m, I, J);
    std::vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 42);
    constraintnp= MSSP_unweighted_framework(g1, closest_pivots,2,links,sp);
    for(int i=0;i<closest_pivots.size();i++){
        pivots[closest_pivots[i]]=1;
    }

    std::set<int> pset;
    for(int i=0;i<closest_pivots.size();i++){
        pset.insert(closest_pivots[i]);
    }
    std::vector<int> pivot;
    /*randperm(p,n,pivot);
    sort(pivot.begin(),pivot.end());*/
    for(std::set<int>::iterator it=pset.begin();it!=pset.end();it++){
        pivot.push_back(*it);
        // std::cout<<*it<<" ";
    }
    for(int i=0;i<pivot.size();i++){
        for(int j=0;j<n;j++){
            if(pivot[i]!=j){
                float w = 1.0f / (sp[pivot[i]][j] * sp[pivot[i]][j]);
                constraints.push_back(constraint(pivot[i],j,sp[pivot[i]][j],w,pa));
            }
        }
    }
    /*for(int i=0;i<constraintnp.size();i++){
        consSGD.insert(constraintnp[i]);
    }*/

}
void graph::append_Np_range(force pa) {

    //consSGD.insert()
    std::vector<constraint> constraintnp;
    int p=100;
    int I[m];
    int J[m];
    for (int i = 0; i < m; i++) {
        I[i]=links[i].source;
        J[i]=links[i].target;
    }
    pivots.resize(n);
    for(int i=0;i<pivots.size();i++){
        pivots[i]=0;
    }
    std::vector<std::vector<int> > g1 = build_graph_unweighted(n, m, I, J);
    std::vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 42);
    constraintnp= MSSP_unweighted_framework(g1, closest_pivots,2,links,sp);
    for(int i=0;i<closest_pivots.size();i++){
        pivots[closest_pivots[i]]=1;
    }

    std::set<int> pset;
    for(int i=0;i<closest_pivots.size();i++){
        pset.insert(closest_pivots[i]);
    }
    std::vector<int> pivot;
    /*randperm(p,n,pivot);
    sort(pivot.begin(),pivot.end());*/
    for(std::set<int>::iterator it=pset.begin();it!=pset.end();it++){
        pivot.push_back(*it);
        // std::cout<<*it<<" ";
    }
    for(int i=0;i<pivot.size();i++){
        for(int j=0;j<n;j++){
            if(pivot[i]!=j){
                float w = 1.0f / (sp[pivot[i]][j] * sp[pivot[i]][j]);
                //constraints.push_back(constraint(pivot[i],j,sp[pivot[i]][j],w,pa));
                int src,tgt;
                if(pivot[i]>j){
                    tgt=pivot[i],src=j;
                }else{
                    src=pivot[i],tgt=j;
                }
                if(append_cons[src][tgt].flag){
                append_cons[src][tgt].force_m.push_back(pa);
                }else{
                append_cons[src][tgt].flag=true;
                append_cons[src][tgt].d=sp[src][tgt];
                append_cons[src][tgt].w=w;
                append_cons[src][tgt].force_m.push_back(pa);
                }
                using_NP=true;
                // auto it = std::set_sgd.find(constraint_sgd(src,tgt));
                // if(it!=std::set_sgd.end()){
                //     it->force_sgd.push_back(pa);
                // }else{
                //     std::vector<force> f;
                //     f.push_back(pa);
                //     std::set_sgd.insert(constraint_sgd(src,tgt,sp[src][tgt],w,f));
                // }
                
            }
        }
    }
    /*for(int i=0;i<constraintnp.size();i++){
        consSGD.insert(constraintnp[i]);
    }*/

}
void graph::append_S_range(int neighbor,para pa) {
    //consSGD.insert()
    std::vector<Edge> e[n];
    int st, end, weight;
    for (int i=0; i < m; i++)
    {
        st=links[i].source;
        end=links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);    //
        st=links[i].target;
        end=links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    for(int i=0;i<n;i++){
        std::vector<int> kn;
        kn.push_back(i);
        DFS2(n,m,e,i,kn,neighbor);
        for(int j=0;j<kn.size();j++){
            if(kn[j]<i){
                float w = 1.0f / (sp[i][kn[j]] * sp[i][kn[j]]);
                constraints.push_back(constraint(i,kn[j],sp[i][kn[j]],w,pa));
            }
        }
    }
}
void graph::append_S_range(int neighbor,force pa) {
    //consSGD.insert()
    std::vector<Edge> e[n];
    int st, end, weight;
    for (int i=0; i < m; i++)
    {
        st=links[i].source;
        end=links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
        st=links[i].target;
        end=links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    for(int i=0;i<n;i++){
        std::vector<int> kn;
        kn.push_back(i);
        DFS2(n,m,e,i,kn,neighbor);
        for(int j=0;j<kn.size();j++){
            if(kn[j]<i){
                float w = 1.0f / (sp[i][kn[j]] * sp[i][kn[j]]);
                int src=i,tgt=kn[j];
                 if(append_cons[src][tgt].flag){
                append_cons[src][tgt].force_m.push_back(pa);
                 }else{
                append_cons[src][tgt].flag=true;
                append_cons[src][tgt].d=sp[src][tgt];
                append_cons[src][tgt].w=w;
                append_cons[src][tgt].force_m.push_back(pa);
                }
                // auto it = std::set_sgd.find(constraint_sgd(i,kn[j]));
                // if(it!=std::set_sgd.end()){
                //     it->force_sgd.push_back(pa);
                // }else{
                //     std::vector<force> f;
                //     f.push_back(pa);
                //     std::set_sgd.insert(constraint_sgd(i,kn[j],sp[i][kn[j]],w,f));
                // }
                
            }
        }
    }
}
bool graph::power_law_graph() {
    std::vector<int> mask(m+1,0);
    std::vector<int> count(n,0);
    int maxD=0;
    bool isPow=false;
    for(int i=0;i<m;i++){
        count[links[i].source]++;
        count[links[i].target]++;
    }
    for(int i=0;i<this->n;i++){
        mask[count[i]]++;
        maxD=std::max(maxD,mask[count[i]]);
    }
    if(mask[1]>0.8*maxD&&mask[1]>0.3*m) {isPow= true;std::cout<<"graph is low power graph,std::set p=1.8"<<std::endl;}
    return isPow;
}

void graph::preSolveSGD(double eps,int t_max,int seed,float eta_max){
    float w_min = 100000000;
    float w_max = 0;
    for(int k=0;k<constraints.size();k++){
        float tempw=constraints[k].w;
        w_min = std::min(tempw, w_min);
        w_max = std::max(tempw, w_max);
    }

    if(paraSGD.pmds==2){
        int p=100;
        if(p>n) p=n;
        int I[m];
        int J[m];
        for (int i = 0; i < m; i++) {
            I[i]=links[i].source;
            J[i]=links[i].target;
        }
        std::vector<std::vector<int> > g1 = build_graph_unweighted(n, m, I, J);
        std::vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 48);
        std::set<int> pset;
        for(int i=0;i<closest_pivots.size();i++){
            pset.insert(closest_pivots[i]);
        }
        std::vector<int> pivot;
        //randperm(p,n,pivot);
        //sort(pivot.begin(),pivot.end());
        for(std::set<int>::iterator it=pset.begin();it!=pset.end();it++){
            pivot.push_back(*it);
            // std::cout<<*it<<" ";
        }
        /*for(int i=0;i<pivot.size();i++){
            std::cout<<pivot[i]<<",";
        }*/
        std::cout<<"compute pivotMDS"<<std::endl;
        PivotMDS(nodes,pivot.size(),n,sp,pivot);
        float maxx=0,minx=100000000000,maxy=0,miny=10000000000;

        for(int i=0;i<nodes.size();i++){
            if(nodes[i].x>maxx) maxx = nodes[i].x;
            if(nodes[i].x<minx) minx = nodes[i].x;
            if(nodes[i].y>maxy) maxy = nodes[i].y;
            if(nodes[i].y<miny) miny = nodes[i].y;
        }
        float lenth=std::max(maxx-minx,maxy-miny);
        for(int i=0;i<nodes.size();i++){
            nodes[i].x-=minx;
            nodes[i].y-=miny;
        }
        for(int i=0;i<nodes.size();i++){
            nodes[i].x/=lenth;
            nodes[i].y/=lenth;
        }
    }

    rk_seed(seed, &rstate);

    fisheryates_shuffle(constraints, rstate);
    float epsilon = 0.1;
    //eta_max = 100.0f;//1.0f/w_min;
    float eta_min = eps;//eps/w_max;
    float lambd = log(eta_min / eta_max) / (t_max - 1);
    for (int i = 0; i < t_max; i++) {
        schedule.push_back(eta_max * exp(lambd * i));
    }
    float eta_maxt = 1.0f / w_min;
    float eta_mint = epsilon / w_max;

    float lambdt = log(eta_mint / eta_maxt) / (t_max - 1);
    for (int i = 0; i < t_max; i++) {
        schedulet.push_back(eta_maxt * exp(lambdt * i));
        //std::cout<<"schedule"<<schedule[i]<<std::endl;
    }
}
void graph::preSolveSGD2(double eps,int t_max,int seed,float eta_max){

    //constraints_sgd.assign(std::set_sgd.begin(),std::set_sgd.end());
    int num_of_node=nodes.size();
    if(t_schema==false){
        for(int i=0;i<num_of_node;i++){
            for(int j=0;j<num_of_node;j++){
                if(append_cons[i][j].flag){
                 constraints_sgd.push_back(constraint_sgd(i,j,append_cons[i][j].d,append_cons[i][j].w,append_cons[i][j].force_m));
                 }
            }
    }
    }else{
        fprintf(stderr,"using t force\n");
        for(int i=0;i<num_of_node;i++){
            for(int j=0;j<num_of_node;j++){
                if(append_cons[i][j].flag){
                 constraints_sgd.push_back(constraint_sgd(i,j,append_cons[i][j].d,append_cons[i][j].w,append_cons[i][j].t_force_m));
                 }
            }
        }
    }
    
    std::cout<<"cons size="<<constraints_sgd.size()<<std::endl;

    float w_min = 100000000;
    float w_max = 0;
    for(int k=0;k<constraints_sgd.size();k++){
        float tempw=constraints_sgd[k].w;
        w_min = std::min(tempw, w_min);
        w_max = std::max(tempw, w_max);
    }
    rk_seed(seed, &rstate);
    fisheryates_shuffle(constraints_sgd, rstate);
    //float epsilon = 0.1;
    //eta_max = 100.0f;//1.0f/w_min;
    float eta_min = eps;//eps/w_max;
    float lambd = log(eta_min / eta_max) / (t_max - 1);
    for (int i = 0; i < t_max; i++) {
        schedule.push_back(eta_max * exp(lambd * i));
    }
    // float eta_maxt = 1.0f / w_min;
    // float eta_mint = epsilon / w_max;

    // float lambdt = log(eta_mint / eta_maxt) / (t_max - 1);
    // for (int i = 0; i < t_max; i++) {
    //     schedulet.push_back(eta_maxt * exp(lambdt * i));
    // }
}
 void graph::preSolveSGD_linear(float e_max,float e_min, int t_max , int seed){
    int num_of_node=nodes.size();
    if(t_schema==false){
        for(int i=0;i<num_of_node;i++){
            for(int j=0;j<num_of_node;j++){
                if(append_cons[i][j].flag){
                 constraints_sgd.push_back(constraint_sgd(i,j,append_cons[i][j].d,append_cons[i][j].w,append_cons[i][j].force_m));
                 }
            }
    }
    }else{
        fprintf(stderr,"using t force\n");
        for(int i=0;i<num_of_node;i++){
            for(int j=0;j<num_of_node;j++){
                if(append_cons[i][j].flag){
                 constraints_sgd.push_back(constraint_sgd(i,j,append_cons[i][j].d,append_cons[i][j].w,append_cons[i][j].t_force_m));
                 }
            }
        }
    }
    
    std::cout<<"cons size="<<constraints_sgd.size()<<std::endl;
    // float w_min = 100000000;
    // float w_max = 0;
    // for(int k=0;k<constraints_sgd.size();k++){
    //     float tempw=constraints_sgd[k].w;
    //     w_min = min(tempw, w_min);
    //     w_max = max(tempw, w_max);
    // }
    rk_seed(seed, &rstate);
    fisheryates_shuffle(constraints_sgd, rstate);
    float epsilon = 0.1;
    //eta_max = 100.0f;//1.0f/w_min;
   // eta_min = eps;//eps/w_max;
    //float lambd = log(eta_min / eta_max) / (t_max - 1);
    double alpha=1.0f;

    alphaDecay=1-pow(e_min,1.0/(float)(t_max));
    for (int i = 0; i < t_max; i++) {
        schedule.push_back(alpha);
        alpha*=(1-alphaDecay);
    }
 }

bool graph::solveSGD(int iter) {
    //SGD解法

    double Delta=0.001;
    double delta_avg=0;
    double sum_delta=0;
    float alpha=0;

    if (iter % 10 == 1)fisheryates_shuffle(constraints, rstate);
    //clock_t time_end = clock();
    //std::cout << "time ran:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << " " <<iter<< std::endl;
    double Delta_max=0;
    // std::cout<<"11111"<<std::endl;
//#pragma omp parallel for schedule(dynamic)
//std::cout<<nodes[0].x<<","<<nodes[0].y<<std::endl;
    for (int con = 0; con < constraints.size(); con++) {
        float rx = 0;
        float ry = 0;

        const constraint &t = constraints[con];
        const float &q = t.pa.q;
        const float &p = t.pa.p;
        const float &a = t.pa.a;
        const float &r = t.pa.r;
        const float &ka = t.pa.ka;
        const float &kr = t.pa.kr;
        /*if(iter==0){
            std::cout<<"q="<<q<<"p="<<p<<"a="<<a<<"r="<<r<<"ka="<<ka<<"kr="<<kr<<std::endl;
        }*/
        const float &pm=t.pa.pmds;
        const int &i = t.i, &j = t.j;
        const double &w_ij = t.w;
        const double &d_ij = t.d;
        if(d_ij==0) continue;
        alpha = schedule[iter];
        //std::cout<<con<<" get schedule"<<std::endl;
        //if(alpha>1.0f)alpha=1.0f;
        float mvx = (nodes[i].x - nodes[j].x);
        float mvy = (nodes[i].y - nodes[j].y);
        float dist = sqrtf(mvx * mvx + mvy * mvy);
        //std::cout<<con<<" get dist"<<std::endl;
        if (dist < 0.1f)dist = 0.1f;
        float fa=0;
        float fccur = 1 / dist;
        float wcur = 1 / d_ij;
        if(ka!=0){
            fa=ka;
            //fa=pow(dist,q);
            //fa*=pow(d_ij,a);
            if(q-(int)q==0){
                if(q<0){
                    for (int pi = 0; pi < -q; pi++)
                        fa *= fccur;
                }else if(q>0){
                    for (int pi = 0; pi < q; pi++)
                        fa *= dist;
                }
            }else{
                //std::cout<<"q="<<q<<std::endl;
                fa*=pow(dist,q);
            }
            if(a-(int)a==0){
                if (a > 0) {
                    for (int pi = 0; pi < a; pi++)
                        fa *= d_ij;
                }else if (a < 0) {
                    for (int pi = 0; pi < -a; pi++)
                        fa *= wcur;
                }
            }else{
                //std::cout<<"a="<<a<<std::endl;
                fa*=pow(d_ij,a);
            }
        }
        //a
        float fr=0;
        if(kr!=0){
            fr=kr;
            if(p-(int)p==0){
                if(p<0){
                    for (int pi = 0; pi < -p; pi++)
                        fr *= fccur;
                }else if(p>0){
                    for (int pi = 0; pi < p; pi++)
                        fr *= dist;
                }
            }else{
                //std::cout<<"p="<<p<<std::endl;
                fr*=pow(dist,p);
            }
            //fr=pow(d_ij,r);
            //fr*=pow(dist,p);
            if(r-(int)r==0){
                if (r > 0) {
                    for (int pi = 0; pi < r; pi++)
                        fr *= d_ij;
                }else if (r < 0) {
                    for (int pi = 0; pi < -r; pi++)
                        fr *= wcur;
                }
            }else{fr=pow(d_ij,r);}
            if(fr>d_ij) fr=d_ij;
        }
        //r
        float f = alpha * (fa - fr);
        if(ka==0||kr==0){
           /* if (abs(f) > 0.999 * abs(dist)) {
                //f = dist;
                f = dist*0.99;
            }*/
            if (abs(f) > abs(dist)) {
                //f = dist;
                f = dist;
            }
        }else{
            if (abs(f) > abs(dist-d_ij)) {
                f = dist-d_ij;
            }
        }
        rx += 1.0f * mvx / dist * f * 0.5;
        ry += 1.0f * mvy / dist * f * 0.5;

            //if(iter==0) std::cout<<"no pivot"<<std::endl;
            nodes[i].x -= rx;
            nodes[i].y -= ry;
            nodes[j].x += rx;
            nodes[j].y += ry;

    }
    return 1;
}
bool graph::solveSGD2(int iter) {
    //SGD解法
    double Delta=0.001;
    double delta_avg=0;
    double sum_delta=0;
    float alpha=0;
    if (iter % 10 == 1)fisheryates_shuffle(constraints_sgd, rstate);

    double Delta_max=0;
    // std::cout<<"11111"<<std::endl;
//#pragma omp parallel for schedule(dynamic)
//std::cout<<nodes[0].x<<","<<nodes[0].y<<std::endl;
    for (int con = 0; con < constraints_sgd.size(); con++) {
        float rx = 0;
        float ry = 0;

        const constraint_sgd &t = constraints_sgd[con];
        const int &i = t.i, &j = t.j;
        const double &w_ij = t.w;
        const double &d_ij = t.d;
        alpha = schedule[iter];
        /*if(iter==0){
            std::cout<<"i="<<i<<"j="<<j<<"w="<<w_ij<<"d="<<d_ij<<std::endl;
        }*/

        if(d_ij==0) continue;

        float mvx = (nodes[i].x - nodes[j].x);
        float mvy = (nodes[i].y - nodes[j].y);
        float dist = sqrtf(mvx * mvx + mvy * mvy);

        if (dist < 0.1f)dist = 0.1f;
        float fccur = 1 / dist;
        float wcur = 1 / d_ij;
        float f=0;
        bool isSping=true;
        for(int _k=0;_k<t.force_sgd.size();_k++){
            const float &k = t.force_sgd[_k].k;
            const float &a = t.force_sgd[_k].a;
            const float &b = -1.0f*t.force_sgd[_k].b;
            /* if(con==0&&iter==0)
             std::cout<<"k="<<k<<"a="<<a<<"b="<<b<<std::endl;*/
            float _f=0;
            _f=k;
            if(k<0) isSping= false;
            if(a-(int)a==0){
                if(a<0){
                    for (int pi = 0; pi < -a; pi++)
                        _f *= fccur;
                }else if(a>0){
                    for (int pi = 0; pi < a; pi++)
                        _f *= dist;
                }
            }else{
                _f*=pow(dist,a);
            }
            if(b-(int)b==0){
                if (b > 0) {
                    for (int pi = 0; pi < b; pi++)
                        _f *= d_ij;
                }else if (b < 0) {
                    for (int pi = 0; pi < -b; pi++)
                        _f *= wcur;
                }
            }else{
                //std::cout<<"a="<<a<<std::endl;
                _f*=pow(d_ij,b);
            }
            //if(_f<-1*d_ij){_f=-1*d_ij;}
            f+=_f;
        }
        f*=alpha;
        if(isSping){
            if (abs(f) > 0.999 * abs(dist)) {
                f = dist*0.99;
            }
        }else{
            if (abs(f) > abs(dist-d_ij)) {
                f = dist-d_ij;
            }
        }
       // std::cout<<"isSpring="<<isSping<<std::endl;
        rx += 1.0f * mvx / dist * f * 0.5;
        ry += 1.0f * mvy / dist * f * 0.5;

            nodes[i].x -= rx;
            nodes[i].y -= ry;
            nodes[j].x += rx;
            nodes[j].y += ry;

    }
    return 1;
}
bool graph::solveSGD_t(int iter) {
    //SGD解法
    double Delta=0.001;
    double delta_avg=0;
    double sum_delta=0;
    float alpha=0;
    if (iter % 10 == 1)fisheryates_shuffle(constraints_sgd, rstate);

    double Delta_max=0;
    for (int con = 0; con < constraints_sgd.size(); con++) {
        float rx = 0;
        float ry = 0;

        const constraint_sgd &t = constraints_sgd[con];
        const int &i = t.i, &j = t.j;
        const double &w_ij = t.w;
        const double &d_ij = t.d;
        alpha = schedule[iter];
        /*if(iter==0){
            std::cout<<"i="<<i<<"j="<<j<<"w="<<w_ij<<"d="<<d_ij<<std::endl;
        }*/

        if(d_ij==0) continue;

        float mvx = (nodes[i].x - nodes[j].x);
        float mvy = (nodes[i].y - nodes[j].y);
        float dist = sqrtf(mvx * mvx + mvy * mvy);
        float eps=0.1f;
        dist=dist+eps;
       // if (dist < 0.1f)dist = 0.1f;
        float fccur = 1 / dist;
        float wcur = 1 / d_ij;
        float f=0;
        float t_dist=1/(1+dist*dist);
        bool isSping=true;
        for(int _k=0;_k<t.force_t.size();_k++){
            const float &k = t.force_t[_k].k;
            const float &a = t.force_t[_k].a;
            const float &b = -1.0f*t.force_t[_k].b;
            const float &_t=t.force_t[_k].t;
        
            float _f=0;
            _f=k;
            if(k<0) isSping= false;
            if(a-(int)a==0){
                if(a<0){
                    for (int pi = 0; pi < -a; pi++)
                        _f *= fccur;
                }else if(a>0){
                    for (int pi = 0; pi < a; pi++)
                        _f *= dist;
                }
            }else{
                _f*=pow(dist,a);
            }
            if(b-(int)b==0){
                if (b > 0) {
                    for (int pi = 0; pi < b; pi++)
                        _f *= d_ij;
                }else if (b < 0) {
                    for (int pi = 0; pi < -b; pi++)
                        _f *= wcur;
                }
            }else{
               
                _f*=pow(d_ij,b);
            }
            if(_t-(int)_t==0){
                if (_t > 0) {
                    for (int pi = 0; pi < _t; pi++){
                        _f *= t_dist; 
                    }
                        
                }else if (_t < 0) {
                    for (int pi = 0; pi < -_t; pi++)
                        _f *= (1.0/t_dist);
                }
            }else{
                _f*=pow(t_dist,_t);
            }
            //if(_f<-1*d_ij){_f=-1*d_ij;}
            f+=_f;
        }
        f*=alpha;
        //if(isSping){
            if (abs(f) > 0.999 * abs(dist)) {
                if(f>0){
                    f = dist*0.999;  
                }else{
                    f = -1.0*dist*0.999;
                }
                
             }
        // }
       // else{
            // if (abs(f) > abs(dist-d_ij)) {
            //     if(f>0){
            //         f = dist-d_ij;
            //     }else{
            //         f=d_ij-dist;
            //     }
                
            // }
       // }
        rx += 1.0f * mvx / dist * f * 0.5;
        ry += 1.0f * mvy / dist * f * 0.5;

            nodes[i].x -= rx;
            nodes[i].y -= ry;
            nodes[j].x += rx;
            nodes[j].y += ry;

    }
    return 1;
}

void graph::append_BH_range(int neighbor,para pa) {
    //consSGD.insert()
    std::vector<Edge> e[n];
    int st, end, weight;
    for (int i=0; i < m; i++)
    {
        st=links[i].source;
        end=links[i].target;
        Edge tmp;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
        st=links[i].target;
        end=links[i].source;
        tmp.w = DEFAULT_WEIGHT;
        tmp.end = end;
        e[st].push_back(tmp);
    }
    for(int i=0;i<n;i++){
        std::vector<int> kn;
        kn.push_back(i);
        DFS2(n,m,e,i,kn,neighbor);
        for(int j=0;j<kn.size();j++){
            if(kn[j]<i&&sp[i][kn[j]]!=0){
                float w = 1.0f / (sp[i][kn[j]] * sp[i][kn[j]]);
                constraintbh.push_back(constraint(i,kn[j],sp[i][kn[j]],w,pa));
            }
        }
    }
}
void graph::append_BH_range2(force pa) {
    bh_forces.push_back(pa);
}
void graph::append_BH_t_force(t_force pa){
    bh_t_forces.push_back(pa);
}
void graph::preSolveBH(){
    for(int i=0;i<n;i++){
        count.push_back(0);
        force_r.push_back(0);
        force_r.push_back(0);
    }
    for(int k=0;k<constraintbh.size();k++){
        count[constraintbh[k].i]++;
        count[constraintbh[k].j]++;
    }
}
void graph::preSolveBH2(){
    for(int i=0;i<n;i++){
        force_r.push_back(0);
        force_r.push_back(0);
    }
}
void graph::preSolve_t_BH(int t_max,double decay,double _threshold){
    for(int i=0;i<n;i++){
        force_r.push_back(0);
        force_r.push_back(0);
    }
   //
    alphaDecay=1-pow(_threshold,1.0/(float)(t_max));
    //alphaDecay=decay;
    threshold=_threshold;
}
 void graph::solve_t_BH(int iter){
    
    for(int _i=0;_i<bh_t_forces.size();_i++){
        const float &kr=bh_t_forces[_i].k;
        const float &p=bh_t_forces[_i].a;
        const float &tf=bh_t_forces[_i].t;
        float alpha=1.0f;
        float C = kr, K = 1, linkLen = 1,theta = 1.0;
        //斥力做第三个力用0.01 斥力做第二个力用1
       // std::cout<<"alphadecay="<<alphaDecay<<std::endl;
        alpha=pow((1-alphaDecay),iter);
        if(alpha<threshold) alpha=threshold;
        //std::cout<<force[0]<<" "<<force[1]<<" ";
        Quad* QuadTree = new Quad(nodes, n);
        int nn=n;
        float dp=p;
   // #pragma omp parallel for schedule(guided)
        for (int ii = 0; ii < nn; ii++) {
            Quad_t_Force(QuadTree, force_r[ii * 2], force_r[ii * 2 + 1], nodes[ii].x, nodes[ii].y, theta, alpha * C, ii,dp,tf);
        }
        delete QuadTree;
        for (int ii = 0; ii < n; ii++) {
            nodes[ii].x += force_r[2 * ii]*= velocityDecay ;
            nodes[ii].y += force_r[2 * ii + 1]*= velocityDecay;//  
        }
    }
}
void graph::solveBH(int iter) {
    const float &q = paraBH.q;
    const float &p = paraBH.p;
    const float &a = paraBH.a;
    const float &r = paraBH.r;
    const float &ka = paraBH.ka;
    const float &kr = paraBH.kr;
    //BH解法
    float alpha=1.0;
    float C = -1*kr, K = 1, linkLen = 1,theta = 1.0;
    //斥力做第三个力用0.01 斥力做第二个力用1
    alpha=pow((1-alphaDecay),iter);
    if(alpha<0.0001f) alpha=0.0001f;
    //std::cout<<force[0]<<" "<<force[1]<<" ";
    Quad* QuadTree = new Quad(nodes, n);
    int nn=n;
    float dp=p;
//#pragma omp parallel for schedule(guided)
    for (int ii = 0; ii < nn; ii++) {
        QuadForce(QuadTree, force_r[ii * 2], force_r[ii * 2 + 1], nodes[ii].x, nodes[ii].y, theta, alpha * C, ii,dp);
    }

    delete QuadTree;
        for (int ii = 0; ii < n; ii++) {
            nodes[ii].x += force_r[2 * ii]*= velocityDecay ;
            nodes[ii].y += force_r[2 * ii + 1] *= velocityDecay ;//*= velocityDecay
        }
}
void graph::solveBH2(int iter) {
    for(int _i=0;_i<bh_forces.size();_i++){
        const float &kr=bh_forces[_i].k;
        const float &p=-1.0*bh_forces[_i].a;
        float alpha=1.0;
        float C = kr, K = 1, linkLen = 1,theta = 1.0;
        //斥力做第三个力用0.01 斥力做第二个力用1
        alpha=pow((1-alphaDecay),iter);
        if(alpha<0.0001f) alpha=0.0001f;
        //std::cout<<force[0]<<" "<<force[1]<<" ";
        Quad* QuadTree = new Quad(nodes, n);
        int nn=n;
        float dp=p;
//#pragma omp parallel for schedule(guided)
        for (int ii = 0; ii < nn; ii++) {
            QuadForce(QuadTree, force_r[ii * 2], force_r[ii * 2 + 1], nodes[ii].x, nodes[ii].y, theta, alpha * C, ii,dp);
        }
        delete QuadTree;
        for (int ii = 0; ii < n; ii++) {
            nodes[ii].x += force_r[2 * ii]*= velocityDecay ;
            nodes[ii].y += force_r[2 * ii + 1] *= velocityDecay ;
        }
    }

}



/*
void graph::outToFile(std::string filename,std::string methodname) {

    std::string str1= ".\\" + methodname;
    //system(str1.c_str());
    if (access(str1.c_str(),6)==-1)
    {
        mkdir(str1.c_str());
    }
    std::string str2= ".\\" + methodname+"\\"+ filename.substr(0, filename.length() - 4);
    if (access(str2.c_str(),6)==-1)
    {
        mkdir(str2.c_str());
    }
    std::string filename1 = methodname+"\\"+ filename.substr(0, filename.length() - 4)+"\\"+"nodes.csv";
    freopen((char*) filename1.data(),"w",stdout);
    //freopen(str,"w",stdout);
    std::cout<<"positionx"<<","<<"positiony"<<std::endl;

    for(int i=0;i<this->n;i++){
        std::cout<<nodes[i].x<<","<<nodes[i].y<<std::endl;

    }

    //freopen("/dev/tty", "w", stdout);
    freopen("CON", "r", stdin);
    freopen("CON", "w", stdout);
}

void graph::RelationOutToFile(std::string filename,std::string methodname) {

    std::string str1= ".\\" + methodname;
    //system(str1.c_str());
    if (access(str1.c_str(),6)==-1)
    {
        mkdir(str1.c_str());
    }
    std::string str2= ".\\" + methodname+"\\"+ filename.substr(0, filename.length() - 4);
    if (access(str2.c_str(),6)==-1)
    {
        mkdir(str2.c_str());
    }
    std::string filename1 = methodname+"\\"+ filename.substr(0, filename.length() - 4)+"\\"+"links.csv";
    freopen((char*) filename1.data(),"w",stdout);
    std::cout<<"source"<<","<<"target"<<std::endl;
    for(int i=0;i<this->m;i++){
        std::cout<<links[i].source<<","<<links[i].target<<std::endl;
    }
    //freopen("/dev/tty", "w", stdout);
    freopen("CON", "r", stdin);
    freopen("CON", "w", stdout);
}
*/
void graph::initPosition(std::string filename) {
    std::ifstream input(filename);
    std::string s;
    //std::cout<<"read data"<<std::endl;
    if(nodes.size()!=n){
        nodes.resize(n);
    }
    int index=0;

    while(getline(input, s)){
        std::string s1,s2;
        std::istringstream is(s);
        is>>s1>>s2;
        nodes[index].x= str2float(s1);
        nodes[index].y= str2float(s2);
        index++;
    }
}
void graph::drawSVG(std::string filename, std::string method,bool output_pos) {
    std::string filename1 = method+"_" +filename.substr(0, filename.length() - 4)+"_layout.svg";
    float maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9;
    for(int i=0;i<nodes.size();i++){
        if(nodes[i].x>maxx) maxx = nodes[i].x;
        if(nodes[i].x<minx) minx = nodes[i].x;
        if(nodes[i].y>maxy) maxy = nodes[i].y;
        if(nodes[i].y<miny) miny = nodes[i].y;
    }
    float lenth=std::max(maxx-minx,maxy-miny);
    for(int i=0;i<nodes.size();i++){
        nodes[i].x-=minx;
        nodes[i].y-=miny;
    }
    for(int i=0;i<nodes.size();i++){
        nodes[i].x/=lenth;
        nodes[i].x*=1000;
        nodes[i].x+=100;
        nodes[i].y/=lenth;
        nodes[i].y*=1000+100;
        nodes[i].y+=100;
    }
    fprintf(stderr,"new drawSVG\n");
    // std::vector<int> fflag={0,1,2,3,10,12,15,18,19,21};
    // std::vector<int> f_ind(nodes.size(),0);
    // for(int i=0;i<fflag.size();i++){
    //     f_ind[fflag[i]]=1;
    // }
    //std::vector<int> fflag={72,77,100,101,103,157,257,258,259,260,261,262,263,264,265,267,268,269,271,272,273,274,275,276,277,278,280,281,282,283,284,286,287,288,289,290,291,292,293,294,295,297,298,299,300,301};
    //{265,282,267,72,269,270,264,274,291,159,286,294,303,285,261,301,271,273,284,288,292,100,302,103,280,299,275,263,266,259,283,290,300,77,276,260,293,101,279,296,289,278,262,297,277,295};
    FILE *fp1=freopen((char*) filename1.data(),"w",stdout);
    std::cout<<"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n"
          "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
          "\n"
          "<svg width=\"1200\" height=\"1200\" version=\"1.1\"\n"
          "xmlns=\"http://www.w3.org/2000/svg\">"<<std::endl;
    
    for(int i=0;i<this->m;i++){
        //std::cout<<links[i].source<<","<<links[i].target<<std::endl;
      //  if(f_ind[links[i].source]==0||f_ind[links[i].target]==0){
            std::cout<<"<line x1=\""<<nodes[links[i].source].x<<"\" y1=\""<<nodes[links[i].source].y<<"\" x2=\""<<nodes[links[i].target].x<<"\" y2=\""<<nodes[links[i].target].y<<"\"\n"
             "style=\"stroke:rgb(99,99,99);stroke-width:2\"/>"<<std::endl;
       // }
       
    }
    //  for(int i=0;i<this->m;i++){
    //     //std::cout<<links[i].source<<","<<links[i].target<<std::endl;
    //     if(f_ind[links[i].source]==0||f_ind[links[i].target]==0){
            
    //     }else{
    //          std::cout<<"<line x1=\""<<nodes[links[i].source].x<<"\" y1=\""<<nodes[links[i].source].y<<"\" x2=\""<<nodes[links[i].target].x<<"\" y2=\""<<nodes[links[i].target].y<<"\"\n"
    //          "style=\"stroke:rgb(108,255,219);stroke-width:3\"/>"<<std::endl;
    //     }
       
    // }
    for(int i=0;i<this->n;i++){
        //std::cout<<links[i].source<<","<<links[i].target<<std::endl;
        //if(f_ind[i]==0){
        std::cout<<"<circle cx=\""<<nodes[i].x<<"\" cy=\""<<nodes[i].y<<"\" r=\"3\" stroke=\"green\"\n"
              "stroke-width=\"1.5\" fill=\"blue\"/>"<<std::endl;
        // }else{
        //     std::cout<<"<circle cx=\""<<nodes[i].x<<"\" cy=\""<<nodes[i].y<<"\" r=\"4\" stroke=\"green\"\n"
        //       "stroke-width=\"2.5\" fill=\"red\"/>"<<std::endl;
        // }
    }
    //  for(int i=0;i<fflag.size();i++){
    //     //std::cout<<links[i].source<<","<<links[i].target<<std::endl;
    //     std::cout<<"<circle cx=\""<<nodes[fflag[i]].x<<"\" cy=\""<<nodes[fflag[i]].y<<"\" r=\"5\" stroke=\"black\"\n"
    //           "stroke-width=\"2.5\" fill=\"red\"/>"<<std::endl;
    // }
    std::cout<<"</svg>"<<std::endl;
    if(output_pos){
        std::string filename2 = method+"_" +filename.substr(0, filename.length() - 4)+".csv";
        FILE *fp3=freopen((char*)filename2.data(), "w", stdout);
        std::cout<<"pos_x,pos_y"<<std::endl;
        for(int i=0;i<this->n;i++){
            std::cout<<nodes[i].x<<","<<nodes[i].y<<std::endl;
        }
         std::string filename3 = "link_" +filename.substr(0, filename.length() - 4)+".csv";
        FILE *fp4=freopen((char*)filename3.data(), "w", stdout);
        for(int i=0;i<this->m;i++){
            std::cout<<links[i].source<<","<<links[i].target<<std::endl;
        }
        std::string filename4 =filename.substr(0, filename.length() - 4)+".js";
        FILE *fp5=freopen((char*)filename4.data(), "w", stdout);
        std::cout<<"export function dwt_1005(){"<<std::endl;
        std::cout<<"const s = ["<<std::endl;
        for(int i=0;i<this->m;i++){
            std::cout<<"["<<links[i].source<<","<<links[i].target<<"],"<<std::endl;
        }
        std::cout<<"]"<<std::endl<<"return s;"<<std::endl<<"}"<<std::endl;
    }
     
    //freopen("CON", "r", stdin);
    //freopen("CON", "w", stdout);
     FILE *fp2=freopen("/dev/tty", "w", stdout);
}