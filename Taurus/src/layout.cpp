#include "layout.h"
void SMLayout(graph &g){
    g.initPosition();
    g.solveDij();
    int t_max = 30;
    para pa_sm(1, 1, -2, 1, -1, 0, 1);
    g.paraSGD = pa_sm;
    g.append_N_range(pa_sm);

    g.preSolveSGD(0.01, t_max, 42);
    cout<<"cons size="<<g.constraints.size()<<endl;
    for (int iter = 0; iter < t_max; iter++) {
        bool flag = g.solveSGD(iter);
    }
}
void BSMLayout(graph &g){
    g.initPosition();
    g.solveDij();
    int t_max = 30;
    para pa_new7(1, 1, -1, 1, 1, -1, 1);
    g.paraSGD = pa_new7;

    g.append_N_range(pa_new7);
    g.preSolveSGD(0.01, t_max, 42);
    for (int iter = 0; iter < t_max; iter++) {
        bool flag = g.solveSGD(iter);
    }
}
void FDPLayout(graph &g){
    g.initPosition();
    g.solveDij();
    para pa_fdp2(40, 2, 0, 0, 0, 0, 1);
    para pa_bh(0, 0, 0, 1, 0, 1, 1);// 1/d的斥力
    int t_max=200;
    int c=0;
    t_max = 200;
    g.paraSGD = pa_fdp2;
    g.paraBH = pa_bh;
    g.alphaDecay=0.71;
    c=1;

    g.append_E_range(pa_fdp2);
    g.append_BH_range(0,pa_bh);

    g.preSolveSGD(0.01, t_max, 42);
    g.preSolveBH();
    for (int iter = 0; iter < t_max; iter++) {
        bool flag = g.solveSGD(iter);
        // if(flag==0){cout<<iter<<"break"<<endl;break;}
        if(iter!=(t_max-c)){
            g.solveBH(iter);
        }
    }
}
void LinLogLayout(graph &g){
    g.initPosition();
    g.solveDij();
    para pa_linlog(1, 0, 0, 0, 0, 0, 1);//引力部分
    para pa_linlog2(0, 0, 0, 0.1, 0, 1, 1);//四叉树的参数 kr斥力系数 p斥力的指数 1/d^p p取正
    int t_max = 200;
    g.paraSGD = pa_linlog;
    g.paraBH = pa_linlog2;
    //g.paraBH=pa_bh;
    g.alphaDecay=0.034;
    int c=1;
    g.append_E_range(pa_linlog);
    g.append_BH_range(0,pa_linlog2);
    g.preSolveSGD(0.01, t_max, 42);
    g.preSolveBH();
    for (int iter = 0; iter < t_max; iter++) {
        bool flag = g.solveSGD(iter);
        // if(flag==0){cout<<iter<<"break"<<endl;break;}
        if(iter<(t_max-c)){
            g.solveBH(iter);
        }
    }
}
void FA2Layout(graph &g){
    g.initPosition();
    g.solveDij();
    para pa_fa2(1, 1, 0, 0, 0, 0, 1);
    para pa_bh(0, 0, 0, 1, 0, 1, 1);// 1/d的斥力
    int t_max = 200;
    g.paraSGD = pa_fa2;
    g.paraBH = pa_bh;
    g.deg_rep=true;
    int c=0;
    g.append_E_range(pa_fa2);
    g.append_BH_range(0,pa_bh);
    g.alphaDecay=0.034;

    g.preSolveSGD(0.01, t_max, 42);
    g.preSolveBH();
    for (int iter = 0; iter < t_max; iter++) {
        bool flag = g.solveSGD(iter);
        // if(flag==0){cout<<iter<<"break"<<endl;break;}
        if(iter!=(t_max-c)){
            g.solveBH(iter);
        }
    }
}
void MaxentLayout(graph &g){
    g.initPosition();
    g.solveDij();
    para pa_maxent(1, 1, -2, 1, -1, 0, 2);
    para pa_maxent2(0, 0, 0, 0.01, 0, 1, 1);
    para pa_maxent3(0, 0, 0, 0.0016, 0, 1.8, 1);//for tree
    int t_max = 200;
    g.deg_rep= false;
    g.paraSGD = pa_maxent;
    if(!g.power_law_graph()){
        t_max=200;
        g.paraBH = pa_maxent2;
        g.alphaDecay=0.1;
    }else{
        t_max=200;
        g.paraBH = pa_maxent3;
        g.alphaDecay=0.001;
    }
    g.append_E_range(pa_maxent);
    g.append_BH_range(0,g.paraBH);
    g.preSolveSGD(0.01, t_max, 42,1.0f);
    g.preSolveBH();

    for (int iter = 0; iter < t_max; iter++) {
        bool flag = g.solveSGD(iter);
        g.solveBH(iter);
    }
}
vector<string> Split(const string &s, const string &seperator){
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while(i != s.size()){
        //找到字符串中首个不等于分隔符的字母；
        int flag = 0;
        while(i != s.size() && flag == 0){
            flag = 1;
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[i] == seperator[x]){
                    ++i;
                    flag = 0;
                    break;
                }
        }

        //找到又一个分隔符，将两个分隔符之间的字符串取出；
        flag = 0;
        string_size j = i;
        while(j != s.size() && flag == 0){
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[j] == seperator[x]){
                    flag = 1;
                    break;
                }
            if(flag == 0)
                ++j;
        }
        if(i != j){
            result.push_back(s.substr(i, j-i));
            i = j;
        }
    }
    return result;
}

void Layout(graph &g,vector<force> f){
    //if(num!=f.size()){cout<<"paremeter error!"; return;}
    int num=f.size();
    g.solveDij();
    if(num<2) {fprintf(stderr,"parameter is not enough!\n"); return;}
    cout<<"input parameter is"<<endl;
    for(int i=0;i<f.size();i++){
        cout<<f[i].range<<" "<<f[i].k<<" "<<f[i].a<<" "<<f[i].b<<endl;
    }
    int num_of_node=g.nodes.size();
   /* para pa1(0,0,0,0,0,0,1);
    g.paraSGD=pa1;*/
    bool useQT=false;
    //int count_fr1=0,count_fr2=0;
    bool big_fa=false,count_fr1=false,count_fr2=false;
    for(int i=0;i<num;i++){
        // if(paras[4*i]!=0&&paras[4*i+1]<0){cout<<paras[4*i]<<","<<paras[4*i+1]<<"useQT=false";useQT=false;}
        if(f[i].range==0&&f[i].k<0&&f[i].b!=0){count_fr1=true;}
        if(f[i].range==0&&f[i].k<0&&f[i].b==0){count_fr2=true;}
        if(f[i].range!=1&&f[i].k>0){big_fa= true;}
    }
    //if(big_fa&&count_fr1&&count_fr2){useQT= true;}else if(!big_fa&&!count_fr1)
    if(!big_fa&&!count_fr1&&count_fr2){useQT=true;}
    cout<<"useQT"<<useQT<<endl;
    //malloc
    for(int _i=0;_i<num_of_node;_i++){
        vector<merge_cons> temp_cons(num_of_node);
        g.append_cons.push_back(temp_cons);
    }
    cout<<"append_cons_size="<<g.append_cons.size()<<endl;

    for(int i=0;i<f.size();i++){
        if(!useQT){
            if(f[i].range==0){
                g.append_N_range(f[i]);
                cout<<"append N"<<i<<endl;
            }else if(f[i].range==1){
                g.append_E_range(f[i]);
            }else if(f[i].range==2){
                g.append_S_range(2,f[i]);
            }else{
                g.append_Np_range(f[i]);
            }
        }else if(f[i].range==0&&f[i].k<0){
            g.append_BH_range2(f[i]);
        }else{
            if(f[i].range==0){
                g.append_N_range(f[i]);
            }else if(f[i].range==1){
                g.append_E_range(f[i]);
            }else if(f[i].range==2){
                g.append_S_range(2,f[i]);
            }else{
                g.append_Np_range(f[i]);
            }
        }
    }

    

    int t_max=30;
    //g.alphaDecay=0.71;
    if(useQT){
        g.preSolveBH2();
        t_max=200;
    }
    g.preSolveSGD2(0.01, t_max, 42);
    // cout<<"use"<<useQT<<endl;
    //
    //cout<<"tmax="<<t_max<<endl;
     /*cout<<"cons.size="<<g.constraints_sgd.size()<<"set.size()="<<g.set_sgd.size()<<endl;
    for(int i=0;i<g.constraints_sgd.size();i++){
        cout<<"i="<<g.constraints_sgd[i].i<<",j="<<g.constraints_sgd[i].j<<" ";
    }*/
    for (int iter = 0; iter < t_max; iter++) {
        g.solveSGD2(iter);
        if(useQT){
            g.solveBH2(iter);
        }
    }
}
void MaxentLayout(graph &g,int k_neigh){
    //g.initPosition();
    g.solveDij();
    //para pa_maxent(1, 1, -2, 1, -1, 0, 2);
    force pa_maxent2(0, -0.01, 1, 0);
    force pa_maxent3(0, -0.0016, 1.8, 0);//for tree
    std::vector<force> f;
    f.push_back(force(2,1,1,2));
    f.push_back(force(2, -1,0,1));
    if(!g.power_law_graph()){
        f.push_back(pa_maxent2);
        g.alphaDecay=0.1;
    }else{
        f.push_back(pa_maxent3);
        g.alphaDecay=0.001;
    }

    int t_max = 200;
   
    int num_of_node=g.nodes.size();
    for(int _i=0;_i<num_of_node;_i++){
        vector<merge_cons> temp_cons(num_of_node);
        g.append_cons.push_back(temp_cons);
    }
    cout<<"append_cons_size="<<g.append_cons.size()<<endl;
    g.append_S_range(k_neigh,f[0]);
    //g.append_E_range(f[0]);
    //g.append_E_range(f[1]);
    g.append_S_range(k_neigh,f[1]);
    g.initPivotMDSPosition(200);
    g.append_BH_range2(f[2]);

    g.preSolveBH2();
    g.preSolveSGD2(0.01, t_max, 42,1.0f);
    
    for (int iter = 0; iter < t_max; iter++) {
        g.solveSGD2(iter);
        g.solveBH2(iter);
    }
}
void t_Layout(graph &g,vector<t_force> f){
    //if(num!=f.size()){cout<<"paremeter error!"; return;}
    int num=f.size();
    g.solveDij();
    if(num<2) {fprintf(stderr,"parameter is not enough!\n"); return;}
    cout<<"input parameter is"<<endl;
    for(int i=0;i<f.size();i++){
        cout<<f[i].range<<" "<<f[i].k<<" "<<f[i].a<<" "<<f[i].b<<" "<<f[i].t<<endl;
    }
    int num_of_node=g.nodes.size();
   /* para pa1(0,0,0,0,0,0,1);
    g.paraSGD=pa1;*/
    g.t_schema=true;
    bool useQT=false;
    //int count_fr1=0,count_fr2=0;
    bool big_fa=false,count_fr1=false,count_fr2=false;
    for(int i=0;i<num;i++){
    
        if(f[i].range==0&&f[i].k<0&&f[i].b!=0){count_fr1=true;}
        if(f[i].range==0&&f[i].k<0&&f[i].b==0){count_fr2=true;}
        if(f[i].range!=1&&f[i].k>0){big_fa= true;}
    }
    //if(big_fa&&count_fr1&&count_fr2){useQT= true;}else if(!big_fa&&!count_fr1)
    if(!big_fa&&!count_fr1&&count_fr2){useQT=true;}
    cout<<"useQT"<<useQT<<endl;
    //malloc
    for(int _i=0;_i<num_of_node;_i++){
        vector<merge_cons> temp_cons(num_of_node);
        g.append_cons.push_back(temp_cons);
    }
    cout<<"append_cons_size="<<g.append_cons.size()<<endl;

    for(int i=0;i<f.size();i++){
        if(f[i].range==0){
            if(!useQT){
                 g.append_N_range(f[i]);
                cout<<"append N"<<i<<endl;
            }else{
                g.append_BH_t_force(f[i]);
            }
        }else{
             if(f[i].range==1){
                g.append_E_range(f[i]);}
        }
    }

    int t_max=30;
    //g.alphaDecay=0.71;
    if(useQT){
        t_max=300;
        g.preSolve_t_BH(t_max,0.034,1);
        g.alphaDecay=0.0001;
    }
    g.preSolveSGD2(1, t_max, 42,100);
    //g.preSolveSGD_linear();
 
    cout<<"t_max="<<t_max<<" decay="<<g.alphaDecay<<" threshold="<<g.threshold<<endl;
    for (int iter = 0; iter < t_max; iter++) {
         g.solveSGD_t(iter);
        if(useQT){
            g.solve_t_BH(iter);
        }
       
        
    }
}



