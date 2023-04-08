//
// Created by xml on 2022/8/3.
//
#include <mutex>
#include <cassert>
#include "Multilevel.h"
#ifndef OGDF_MEMORY_POOL_NTS
static std::mutex s_randomMutex;
#endif
static std::mt19937 s_random;
long unsigned int randomSeed()
{
#ifndef OGDF_MEMORY_POOL_NTS
    std::lock_guard<std::mutex> guard(s_randomMutex);
#endif
    return 7*s_random()+3;  // do not directly return seed, add a bit of variation
}

void Multilevel::setSeed(int val)
{
    s_random.seed(val);
}
int Multilevel::randomNumber(int low, int high)
{
    assert(low<=high);
    std::uniform_int_distribution<> dist(low,high);

#ifndef OGDF_MEMORY_POOL_NTS
    std::lock_guard<std::mutex> guard(s_randomMutex);
#endif
    return dist(s_random);
}
void Multilevel::create_multilevel_representations(subgraph g,vector<subgraph> &sg,int &max_level,int seed,int min_Graph_size,int randomTries){
    int bad_edgenr_counter = 0;
    int act_level = 0;
    max_level=30;
    sg[0]=g;
    int st,end,sti,endi;
    for (int i = 0; i < g.numofnode; ++i) {
        vector<adjEdge> tt;
        sg[act_level].edge.push_back(tt);
    }
    sg[act_level].numofnode=g.numofnode;
    sg[act_level].numofedge=g.numofedge;
    for(int i=0;i<g.numofedge;i++){
        st=sg[act_level].sub_edge[i].source;
        end=sg[act_level].sub_edge[i].target;
        adjEdge t_e1,t_e2;
        t_e1.end=end;
        t_e1.length=sg[act_level].sub_edge[i].length;
        t_e2.end=st;
        t_e2.length=sg[act_level].sub_edge[i].length;
        sg[act_level].edge[st].push_back(t_e1);
        sg[act_level].edge[end].push_back(t_e2);
    }
    sg.push_back(g);

    cout<<min_Graph_size<<endl;
    while(sg[act_level].sub_node.size()>min_Graph_size&&edgenumbersum_of_all_levels_is_linear(sg,act_level,bad_edgenr_counter)){
        subgraph a;
        sg.push_back(a);
        init_multilevel_values(sg,act_level);
        partition_galaxy_into_solar_systems(sg,act_level);
        //collaps_solar_systems();
        collaps_solar_systems(sg,act_level);
        act_level++;
        cout<<"act_level============"<<act_level<<endl;
    }
   // cout<<"mm222mm"<<endl;
    max_level=act_level;
   // fprintf(stderr,"act_level=============%d\n",act_level);

}
bool Multilevel::edgenumbersum_of_all_levels_is_linear(vector<subgraph> &sg,int act_level, int &bad_edgenr_counter){
    if(act_level == 0)
        return true;
    else
    {
        if(sg[act_level].sub_edge.size()<=
           0.8 * double (sg[act_level-1].sub_node.size()))
            return true;
        else if(bad_edgenr_counter < 5)
        {
            bad_edgenr_counter++;
            return true;
        }
        return false;
    }
}
void Multilevel::init_multilevel_values(vector<subgraph> &sg,int act_level){
    //
    int num_of_node=sg[act_level].sub_node.size();
    int num_of_edge=sg[act_level].sub_edge.size();
    for (int i = 0; i < num_of_node; ++i) {
        sg[act_level].sub_node[i].type=0;
        sg[act_level].sub_node[i].dedicated_sun_node=-1;
        sg[act_level].sub_node[i].dedicated_sun_distance=0;
        sg[act_level].sub_node[i].dedicated_pm_node=-1;
        sg[act_level].sub_node[i].placed= false;
        sg[act_level].sub_node[i].angle_1=0;
        sg[act_level].sub_node[i].angle_2=2.0*pi;
    }
   /* if(act_level!=0){

    }*/
    for (int i = 0; i < num_of_node;i++){
        int adj_size= sg[act_level].edge[i].size();
        for (int j = 0; j <adj_size ; ++j) {
            int end=sg[act_level].edge[i][j].end;
            sg[act_level].edge[i][j].moon_edge= false;
        }
    }
}
void Multilevel::partition_galaxy_into_solar_systems(vector<subgraph> &sg,int act_level){
    create_suns_and_planets(sg,act_level);
    create_moon_nodes_and_pm_nodes(sg,act_level);
}

void Multilevel::create_suns_and_planets(vector<subgraph> &sg,int act_level){
    vector<int> Node_Set(sg[act_level].sub_node.size(),0);
    vector<int> mass_of_star(sg[act_level].sub_node.size(),0);
    vector<int> planet_nodes;
    vector<int> sun_nodes;
    int node_act_level=sg[act_level].sub_node.size();
    int edge_act_level=sg[act_level].sub_edge.size();
    for(int i=0;i<node_act_level;i++){
        if(act_level==0){
            sg[act_level].sub_node[i].mass=1;
        }
    }
    for(int i=0;i<node_act_level;i++){
       // cout<<"i="<<i<<endl;
        mass_of_star[i]=sg[act_level].sub_node[i].mass;
    }

    for(int i=0;i<edge_act_level;i++){
        int src=sg[act_level].sub_edge[i].source;
        int tgt=sg[act_level].sub_edge[i].target;
       // cout<<"src="<<src<<",tgt"<<tgt<<endl;
        mass_of_star[src]+=sg[act_level].sub_node[tgt].mass;
        mass_of_star[tgt]+=sg[act_level].sub_node[src].mass;
    }
    /*for (int v=0;v<node_act_level;v++){
        fprintf(stderr,"v=%d mass_of_star=%d\n",v,mass_of_star[v]);
    }*/
    //cout<<"aaaaa"<<endl;
    int isNotEmpty=node_act_level;
    int count=0;
    vector<int> S_node;
    vector<int> position_in_node_set;
    for(int i=0;i<node_act_level;i++){
        S_node.push_back(i);
        position_in_node_set.push_back(i);
    }
    int last_selectable_index_of_S_node = node_act_level-1;
    s_random.seed(100);
    while(last_selectable_index_of_S_node>=0){
        vector<int> p_temp;
        planet_nodes.swap(p_temp);
        //cout<<"1111111111111111111111111111"<<endl;
        int rand_index = -1;
        int random_tries=20;
        int temp=get_random_node_with_lowest_star_mass(random_tries,last_selectable_index_of_S_node,mass_of_star,S_node,position_in_node_set);
        //cout<<"sun_node="<<temp<<endl;
        //if(Node_Set[temp]!=0) continue;
        Node_Set[temp]=-1;
        isNotEmpty--;
        sun_nodes.push_back(temp);
        subNode sn;
        //
        sn.id=temp;
        sg[act_level+1].sub_node.push_back(sn);
       // sg[act_level].sub_node[temp].v_higher_level=&sg[act_level+1].sub_node[count];
        sg[act_level+1].sub_node[count].v_lower_level_index=temp;
        sg[act_level].sub_node[temp].v_higher_level_index=count;
        sg[act_level].sub_node[temp].type=1;
        sg[act_level].sub_node[temp].dedicated_sun_node=temp;
        sg[act_level].sub_node[temp].dedicated_sun_distance=0;
        count++;
        for(int i=0;i<sg[act_level].edge[temp].size();i++){
            double dist_to_sun=20.0;
            int tgt=sg[act_level].edge[temp][i].end;
            //if(Node_Set[tgt]!=0)continue;
           // Node_Set[tgt]= -2;
            sg[act_level].sub_node[tgt].type=2;
            sg[act_level].sub_node[tgt].dedicated_sun_node=temp;
           // cout<<"temp="<<temp<<endl;
            sg[act_level].sub_node[tgt].dedicated_sun_distance=dist_to_sun;
            planet_nodes.push_back(tgt);
            //isNotEmpty--;
        }
        for(int i=0;i<planet_nodes.size();i++){
            if(position_in_node_set[planet_nodes[i]] <= last_selectable_index_of_S_node){
                delete_node(planet_nodes[i],last_selectable_index_of_S_node,S_node,position_in_node_set);
            }
        }
        for(int i=0;i<planet_nodes.size();i++){
            for(int j=0;j<sg[act_level].edge[planet_nodes[i]].size();j++){
                int mm=sg[act_level].edge[planet_nodes[i]][j].end;
                if(position_in_node_set[mm] <= last_selectable_index_of_S_node){
                    delete_node(mm,last_selectable_index_of_S_node,S_node,position_in_node_set);
                }
               /* if(Node_Set[mm]==0){
                    Node_Set[mm]= -3;
                    isNotEmpty--;
                }*/
            }
        }
    }
    /*vector<int> sss(node_act_level,0);
    for(int i=0;i<planet_nodes.size();i++){
        cout<<"aaaaaaaaaaaaaa"<<planet_nodes[i]<<endl;
        //sss[sg[act_level].sub_node[i].dedicated_sun_node]=1;
    }*/
   /* int cc=0;
    for (int i = 0; i < node_act_level; ++i) {
        if(sss[i]==1){
            cc+=1;
        }
    }*/

   // cout<<"levelevel"<<endl;
   // cout<<"planet_node"<<endl;
    /*for(int i=0;i<planet_nodes.size();i++){
        cout<<planet_nodes[i]<<",";
    }
    cout<<"sun_node"<<endl;
    for(int i=0;i<sun_nodes.size();i++){
        cout<<sun_nodes[i]<<",";
    }*/
    //init *A_mult_ptr[act_level+1] and set NodeAttributes information for new nodes
    //for(int i=0;i<sg[act_level+1].sub_node.size())
    int num_of_sun=sun_nodes.size();
    //fprintf(stderr,"num_of_sun=%d\n",num_of_sun);
    cout<<"num_of_sun_node="<<num_of_sun<<endl;
    for(int i=0;i<num_of_sun;i++){
        int sun_n=sun_nodes[i];
        //subNode *sn=sg[act_level].sub_node[sun_n].v_higher_level;
        int high_node=sg[act_level].sub_node[sun_n].v_higher_level_index;
        /*sn->weight=sg[act_level].sub_node[sun_n].weight;
        sn->height=sg[act_level].sub_node[sun_n].height;
        sn->m_x=sg[act_level].sub_node[sun_n].m_x;
        sn->m_y=sg[act_level].sub_node[sun_n].m_y;
        sn->v_lower_level=&sg[act_level].sub_node[sun_n];
        sn->v_lower_level_index=sun_n;
        sn->mass=0;*/
        sg[act_level+1].sub_node[high_node].weight=sg[act_level].sub_node[sun_n].weight;
        sg[act_level+1].sub_node[high_node].height=sg[act_level].sub_node[sun_n].height;
        sg[act_level+1].sub_node[high_node].m_x=sg[act_level].sub_node[sun_n].m_x;
        sg[act_level+1].sub_node[high_node].m_y=sg[act_level].sub_node[sun_n].m_y;
       // sg[act_level+1].sub_node[high_node].v_lower_level=&sg[act_level].sub_node[sun_n];
        sg[act_level+1].sub_node[high_node].v_lower_level_index=sun_n;
        sg[act_level+1].sub_node[high_node].mass=0;
    }
}
void Multilevel::delete_node(int del_node,int &last_selectable_index_of_S_node,vector<int> &S_node,vector<int> &position_in_node_set){
    int del_node_index = position_in_node_set[del_node];
    int last_selectable_node = S_node[last_selectable_index_of_S_node];

    S_node[last_selectable_index_of_S_node] = del_node;
    S_node[del_node_index] = last_selectable_node;
    position_in_node_set[del_node] = last_selectable_index_of_S_node;
    position_in_node_set[last_selectable_node] = del_node_index;
    last_selectable_index_of_S_node -=1;
}
int Multilevel::get_random_node_with_lowest_star_mass(int rand_tries,int &last_selectable_index_of_S_node,vector<int> &mass_of_star,vector<int> &S_node,vector<int> &position_in_node_set){
    int rand_index = -1;
    int cmp_mass(0);

    int last_trie_index = last_selectable_index_of_S_node;
    //OGDF_ASSERT(last_trie_index >= 0);

    for (int i = 1; i <= rand_tries && last_trie_index >= 0; ++i) {
        int new_rand_index = randomNumber(0, last_trie_index);
        int mass = mass_of_star[S_node[new_rand_index]];
        get_random_node_common(new_rand_index, last_trie_index,S_node,position_in_node_set);
        if ((i == 1) || (mass<cmp_mass)) {
            rand_index = last_trie_index + 1;
            cmp_mass = mass;
        }
    }
    //OGDF_ASSERT(rand_index != -1);

    return get_random_node_common(rand_index, last_selectable_index_of_S_node,S_node,position_in_node_set);
}
int Multilevel::get_random_node_common(int rand_index, int &last_trie_index,vector<int> &S_node,vector<int> &position_in_node_set)
{
    int random_node = S_node[rand_index];
    int last_trie_node = S_node[last_trie_index];

    S_node[last_trie_index] = random_node;
    S_node[rand_index] = last_trie_node;
    position_in_node_set[random_node] = last_trie_index;
    position_in_node_set[last_trie_node] = rand_index;
    --last_trie_index;
    return random_node;
}
void Multilevel::create_moon_nodes_and_pm_nodes(vector<subgraph> &sg,int act_level){
    int node_act_level=sg[act_level].sub_node.size();
    int edge_act_level=sg[act_level].sub_edge.size();
   cout<<"create_moon_nodes_and_pm_nodes"<<endl;
    for(int i=0;i<node_act_level;i++){
        if(sg[act_level].sub_node[i].type==0){
            //subNode * nearest_neigh_node= nullptr;
            int nearest_neigh_node_index=-1;
            //adjEdge* moon_edge= nullptr;
            double dist_to_nearest_neighbour=0;
            int num_of_adjnode=sg[act_level].edge[i].size();
            int near_node=0;
            for(int j=0;j<num_of_adjnode;j++){
                //cout<<"create_moon_nodes_and_pm_nodes11"<<endl;
                int end=sg[act_level].edge[i][j].end;
                if(i==end) continue;
                 double dist=sg[act_level].edge[i][j].length;
                // cout<<"dist="<<dist<<",";
               int neigh_type= sg[act_level].sub_node[end].type;
              // cout<<"neighbour="<<end<<endl;
                //subNode* neigh_node=&sg[act_level].sub_node[end];
               if((neigh_type==2||neigh_type==3)
               &&(nearest_neigh_node_index== -1||
               dist_to_nearest_neighbour>dist)
               ){
                   // cout<<"type="<<neigh_type<<",";
                   dist_to_nearest_neighbour=dist;
                   nearest_neigh_node_index=end;
                   near_node=j;
               }
            }
           // cout<<endl;
            sg[act_level].edge[i][near_node].moon_edge= true;
            int dedicated_sun_node=sg[act_level].sub_node[nearest_neigh_node_index].dedicated_sun_node;
           // cout<<"near node="<<nearest_neigh_node_index<<" sun node="<<dedicated_sun_node<<",";
            double dedicated_sun_distance = dist_to_nearest_neighbour
                                            + sg[act_level].sub_node[nearest_neigh_node_index].dedicated_sun_distance;
           //cout<<"dedicated_sun_distance"<<dedicated_sun_distance<<" dist_to_nearest_neighbour="<<dist_to_nearest_neighbour<<" nearest_neigh.dedicated_sun_distance="<<sg[act_level].sub_node[nearest_neigh_node_index].dedicated_sun_distance<<endl;
            /*if(nearest_neigh_node_index==7||i==7){
                cout<<"src="<<i<<" tgt="<<nearest_neigh_node_index<<dedicated_sun_distance<<endl;
            }*/
           // cout<<"src="<<i<<" tgt="<<nearest_neigh_node_index<<" distance="<<dedicated_sun_distance<<endl;
            sg[act_level].sub_node[i].type=4;
            sg[act_level].sub_node[i].dedicated_sun_node=dedicated_sun_node;
            sg[act_level].sub_node[i].dedicated_sun_distance=dedicated_sun_distance;
            //sg[act_level].edge[i][j].end
            int end=sg[act_level].edge[i][near_node].end;
            sg[act_level].sub_node[i].dedicated_pm_node=nearest_neigh_node_index;

            sg[act_level].sub_node[nearest_neigh_node_index].type=3;
            sg[act_level].sub_node[nearest_neigh_node_index].moon_List.push_back(i);


        }
    }

}

void Multilevel::collaps_solar_systems(vector<subgraph> &sg,int act_level){
    vector<double> new_edgelength;
   // cout<<"collap1"<<endl;
    calculate_mass_of_collapsed_nodes(sg, new_edgelength,act_level);
   //cout<<"collap2"<<endl;
    create_edges_edgedistances_and_lambda_Lists(sg, act_level);
    //cout<<"collap3"<<endl;
    //delete_parallel_edges_and_update_edgelength(sg, new_edgelength,act_level);
}
void Multilevel::calculate_mass_of_collapsed_nodes(vector<subgraph> &sg,vector<double> &new_edgelength,int act_level){
    int node_act_level=sg[act_level].sub_node.size();
    for(int i=0;i<node_act_level;i++){
        //cout<<"i="<<i<<",";
        int dedicated_sun=sg[act_level].sub_node[i].dedicated_sun_node;
       // cout<<dedicated_sun<<",";
        int high_level_node=sg[act_level].sub_node[dedicated_sun].v_higher_level_index;
        sg[act_level+1].sub_node[high_level_node].mass+=1;
    }
}
void Multilevel::create_edges_edgedistances_and_lambda_Lists(vector<subgraph> &sg,int act_level){
    vector<double> new_edgelength;
    adjEdge e_new;
    vector<adjEdge> inter_solar_system_edges;
    set<adjLink> link_unparallel;
    int edge_act_level=sg[act_level].sub_edge.size();
    for(int i=0;i<edge_act_level;i++){
        int src=sg[act_level].sub_edge[i].source;
        int tgt=sg[act_level].sub_edge[i].target;

        int s_sun_node=sg[act_level].sub_node[src].dedicated_sun_node;
        int t_sun_node=sg[act_level].sub_node[tgt].dedicated_sun_node;

        if(s_sun_node!=t_sun_node){
            /*subNode* high_level_sun_s=sg[act_level].sub_node[s_sun_node].v_higher_level;
            subNode* high_level_sun_t=sg[act_level].sub_node[t_sun_node].v_higher_level;*/
            /*if(s_sun_node>sg[act_level].sub_node.size()||t_sun_node>sg[act_level].sub_node.size()){
                cout<<"out of range     "<<s_sun_node<<" "<<t_sun_node<<",";
            }*/
            int high_level_sun_s=sg[act_level].sub_node[s_sun_node].v_higher_level_index;
            int high_level_sun_t=sg[act_level].sub_node[t_sun_node].v_higher_level_index;
            //cout<<high_level_sun_s<<" "<<high_level_sun_t<<",";
            double length_e=20.0;
            double length_s_edge = sg[act_level].sub_node[src].dedicated_sun_distance;
            double length_t_edge = sg[act_level].sub_node[tgt].dedicated_sun_distance;

            double newlength = length_s_edge + length_e + length_t_edge;
            /*if(high_level_sun_s==7||high_level_sun_t==7){
                cout<<"src="<<high_level_sun_s<<" tgt="<<high_level_sun_t<<" len_s"<<length_s_edge<<" len_t"<<length_t_edge<<endl;
            }*/
          /*  if(high_level_sun_s==7||high_level_sun_t==7){
                fprintf(stderr,"src=%d tgt=%d high_sun_s=%d high_sun_t=%d length_s_edge=%f length_t_edge=%f newlength_e=%f\n",s_sun_node,t_sun_node,high_level_sun_s,high_level_sun_t,length_s_edge,length_t_edge,newlength);
            }*/

            auto it = link_unparallel.find(adjLink(high_level_sun_s,high_level_sun_t));
            auto it_r=link_unparallel.find(adjLink(high_level_sun_t,high_level_sun_s));
            if(it!=link_unparallel.end()||it_r!=link_unparallel.end()){
               // it->force_sgd.push_back(pa);
               //newlength
               if(it!=link_unparallel.end()){
                   it->length+=newlength;
                   it->count++;
               }else{
                   it_r->length+=newlength;
                   it_r->count++;
               }

            }else{
                link_unparallel.insert(adjLink(high_level_sun_s,high_level_sun_t,newlength,1));
            }
            //link_unparallel.insert(adjLink(high_level_sun_s,high_level_sun_t,newlength));

            /*sg[act_level+1].sub_edge.push_back(Link(high_level_sun_s,high_level_sun_t));
            adjEdge ae1,ae2;
            ae1.end=high_level_sun_t;
            ae2.end=high_level_sun_s;
            sg[act_level+1].edge[high_level_sun_s].push_back(ae1);
            sg[act_level+1].edge[high_level_sun_t].push_back(ae2);*/
        }

    }
    sg[act_level+1].sub_edge.assign(link_unparallel.begin(),link_unparallel.end());


    int edge_act_level_1=sg[act_level+1].sub_edge.size();

    for(int i=0;i<sg[act_level+1].sub_node.size();i++){
        vector<adjEdge> tt;
        sg[act_level+1].edge.push_back(tt);
    }
  //  cout<<"create edge"<<endl;
    for(int i=0;i<sg[act_level+1].sub_edge.size();i++){
       // cout<<"edge"<<i<<endl;
        int high_level_sun_s=sg[act_level+1].sub_edge[i].source;
        int high_level_sun_t=sg[act_level+1].sub_edge[i].target;

        if(sg[act_level+1].sub_edge[i].count!=0&&sg[act_level+1].sub_edge[i].count!=1){
            sg[act_level+1].sub_edge[i].length/=sg[act_level+1].sub_edge[i].count;
        }
        //cout<<"edge"<<i<<endl;
        double edge_len=sg[act_level+1].sub_edge[i].length;
       // cout<<"edge"<<i<<endl;
        adjEdge ae1,ae2;
        ae1.end=high_level_sun_t;
        ae1.length=edge_len;
        ae2.end=high_level_sun_s;
        ae2.length=edge_len;
       // cout<<"high_level_sun_t="<<high_level_sun_t<<"high_level_sun_s="<<high_level_sun_s<<endl;
        sg[act_level+1].edge[high_level_sun_s].push_back(ae1);
        sg[act_level+1].edge[high_level_sun_t].push_back(ae2);

    }


    new_edgelength.resize(edge_act_level_1);
    for(int i=0;i<edge_act_level;i++){
        int src=sg[act_level].sub_edge[i].source;
        int tgt=sg[act_level].sub_edge[i].target;
        int s_sun_node=sg[act_level].sub_node[src].dedicated_sun_node;
        int t_sun_node=sg[act_level].sub_node[tgt].dedicated_sun_node;

        double length_e=20.0;
        double length_s_edge = sg[act_level].sub_node[src].dedicated_sun_distance;
        double length_t_edge = sg[act_level].sub_node[tgt].dedicated_sun_distance;
        double newlength = length_s_edge + length_e + length_t_edge;

        double lambda_s = length_s_edge / newlength;
        double lambda_t = length_t_edge / newlength;
        sg[act_level].sub_node[src].lambda.push_back(lambda_s);
        sg[act_level].sub_node[tgt].lambda.push_back(lambda_t);
        sg[act_level].sub_node[src].neighbour_s_node.push_back(t_sun_node);
        sg[act_level].sub_node[tgt].neighbour_s_node.push_back(s_sun_node);

    }
  /*  for(int i=0;i<sg[act_level+1].sub_edge.size();i++){
        int src=sg[act_level+1].sub_edge[i].source;
        int tgt=sg[act_level+1].sub_edge[i].target;
        double len=sg[act_level+1].sub_edge[i].length;
        cout<<" len="<<len<<" src="<<src<<" tgt"<<tgt<<endl;

    }*/

    }

/*void Multilevel::create_initial_placement(vector<subgraph> &sg,int act_level){
    const int BILLION = 1000000000;
    float boxlength=430.0;
    default_random_engine e;
    uniform_int_distribution<int> u(0,BILLION);
    e.seed(time(0));
    int node_of_level=sg[act_level].sub_node.size();
    for(int i=0;i<node_of_level;i++){
        int index=sg[act_level].sub_node[i].id;
        float m_x=double(u(e))/(double)BILLION;
        float m_y=double(u(e))/(double)BILLION;

        sg[0].sub_node[index].m_x=m_x*(boxlength-2)+1;
        sg[0].sub_node[index].m_y=m_y*(boxlength-2)+1;
    }
}*/
void Multilevel::find_initial_placement_for_level(int level,vector<subgraph> &sg){
    vector<int> pm_nodes;
    set_initial_positions_of_sun_nodes(level,sg);
    set_initial_positions_of_planet_and_moon_nodes(level,sg,pm_nodes);
    set_initial_positions_of_pm_nodes(level,sg,pm_nodes);
}
void Multilevel::set_initial_positions_of_sun_nodes(int level,vector<subgraph> &sg){
    int node_num=sg[level+1].sub_node.size();

    for(int i=0;i<node_num;i++){

        int high_node=sg[level+1].sub_node[i].id;
        sg[level].sub_node[sg[level+1].sub_node[i].id].m_x=sg[level+1].sub_node[i].m_x;
        sg[level].sub_node[sg[level+1].sub_node[i].id].m_y=sg[level+1].sub_node[i].m_y;
       // if(level==0){
           // fprintf(stderr,"----level=%d node_id=%d node_mx=%f node_my=%f\n",level,high_node,sg[level+1].sub_node[i].m_x,sg[level+1].sub_node[i].m_y);
       // }
        sg[level].sub_node[sg[level+1].sub_node[i].id].placed= true;
    }
}
void Multilevel::set_initial_positions_of_planet_and_moon_nodes(int level,vector<subgraph> &sg,vector<int> &pm_nodes){
   // float new_pos_x,new_pos_y;
    DPoint new_pos;
    vector<DPoint> L;
    create_all_placement_sectors(level,sg);

    int num_of_node=sg[level].sub_node.size();
    //cout<<"set_p and m nodes"<<endl;
    for (int i = 0; i < num_of_node; ++i) {
        int node_type=sg[level].sub_node[i].type;
        if(node_type==3){
            pm_nodes.push_back(i);
        }else if(node_type==2||node_type==4){
            vector<DPoint> temp;
            L.swap(temp);
            int dedicated_sun=sg[level].sub_node[i].dedicated_sun_node;
            float dedicated_sun_pos_x=sg[level].sub_node[dedicated_sun].m_x;
            float dedicated_sun_pos_y=sg[level].sub_node[dedicated_sun].m_y;
            double dedicated_sun_distance=sg[level].sub_node[i].dedicated_sun_distance;
            DPoint dedicated_sun_pos(dedicated_sun_pos_x,dedicated_sun_pos_y);
            int num_of_edge=sg[level].edge[i].size();
            for (int j = 0; j < num_of_edge; ++j) {
                int end=sg[level].edge[i][j].end;
                if((sg[level].sub_node[i].dedicated_sun_node==
                sg[level].sub_node[end].dedicated_sun_node)&&
                        (sg[level].sub_node[end].type!=1)&&
                        (sg[level].sub_node[end].placed)
                ){
                    new_pos= calculate_position(dedicated_sun_pos,DPoint(sg[level].sub_node[end].m_x,sg[level].sub_node[end].m_y),dedicated_sun_distance,sg[level].edge[i][j].length);
                   /* if(isnan(new_pos.first)){
                        fprintf(stderr,"i=%d j=%d sun=%d dedicated_sun_pos_x=%f dedicated_sun_pos_y=%f sg[level].sub_node[end].m_x=%f sg[level].sub_node[end].m_y=%f distance=%f length=%f\n",i,j,dedicated_sun,dedicated_sun_pos_x,dedicated_sun_pos_y,sg[level].sub_node[end].m_x,sg[level].sub_node[end].m_y,dedicated_sun_distance,sg[level].edge[i][j].length);
                    }*/
                    L.push_back(new_pos);
                }
            }
            if(sg[level].sub_node[i].lambda.empty()){
                if(L.empty()){
                    new_pos= create_random_pos(dedicated_sun_pos,sg[level].sub_node[i].dedicated_sun_distance,sg[level].sub_node[i].angle_1,sg[level].sub_node[i].angle_2);
                    L.push_back(new_pos);

                }
            }else{
                //auto it=sg[level].sub_node[i].lambda_List_ptr->begin();
                int lamda_len=sg[level].sub_node[i].lambda.size();
                int lamda_neigh_len=sg[level].sub_node[i].neighbour_s_node.size();
                for(int k=0;k<lamda_neigh_len;k++){
                    float adj_sun_pos_x,adj_sun_pos_y;
                    int k_neigh=sg[level].sub_node[i].neighbour_s_node.at(k);
                    adj_sun_pos_x=sg[level].sub_node[k_neigh].m_x;
                    adj_sun_pos_y=sg[level].sub_node[k_neigh].m_y;
                    new_pos= get_waggled_inbetween_position(dedicated_sun_pos,DPoint(adj_sun_pos_x,adj_sun_pos_y),sg[level].sub_node[i].lambda.at(k));
                    L.push_back(new_pos);
                }
            }
            sg[level].sub_node[i].placed=true;
            DPoint barycenter_pos= get_barycenter_position(L);
            sg[level].sub_node[i].m_x=barycenter_pos.first;
            sg[level].sub_node[i].m_y=barycenter_pos.second;
           // fprintf(stderr,"planet_and_moon_nodes    m_x=%f m_y=%f\n",sg[level].sub_node[i].m_x,sg[level].sub_node[i].m_y);
        }
    }

}
void Multilevel::set_initial_positions_of_pm_nodes(int level,vector<subgraph> &sg,vector<int> &pm_nodes){
    double moon_dist,lamda;
    int sun_node;
    DPoint sun_pos,moon_pos,new_pos,adj_sun_pos;
    vector<DPoint> L;
    int pm_nodes_len=pm_nodes.size();
    for(int i=0;i<pm_nodes_len;i++){
        vector<DPoint> temp;
        L.swap(temp);
        int v=pm_nodes[i];
        sun_node=sg[level].sub_node[v].dedicated_sun_node;
        sun_pos=DPoint(sg[level].sub_node[sun_node].m_x,sg[level].sub_node[sun_node].m_y);
        double sun_dist=sg[level].sub_node[v].dedicated_sun_distance;
        int adj_len=sg[level].edge[v].size();
        for (int adj  = 0; adj  < adj_len; ++adj ) {
            int end=sg[level].edge[v][adj].end;
            if(end==v)continue;
            if(     (!sg[level].edge[v][adj].moon_edge)&&
                    (sg[level].sub_node[v].dedicated_sun_node==sg[level].sub_node[end].dedicated_sun_node)&&
                    (sg[level].sub_node[end].type!=1)&&
                    (sg[level].sub_node[end].placed)
            ){
                new_pos= calculate_position(sun_pos,DPoint(sg[level].sub_node[end].m_x,sg[level].sub_node[end].m_y),sun_dist,sg[level].edge[v][adj].length);
                //fprintf(stderr,"new_pos_x=%f new_pos_y=%f\n",new_pos.first,new_pos.second);
                L.push_back(new_pos);
            }
        }

        int moon_node_len=sg[level].sub_node[v].moon_List.size();
        for(int moon_node=0;moon_node<moon_node_len;moon_node++){
            int v_moon=sg[level].sub_node[v].moon_List.at(moon_node);
            moon_pos=DPoint (sg[level].sub_node[v_moon].m_x,sg[level].sub_node[v_moon].m_y);
            moon_dist=sg[level].sub_node[v_moon].dedicated_sun_distance;
            lamda=sun_dist/moon_dist;
            new_pos= get_waggled_inbetween_position(sun_pos,moon_pos,lamda);
            L.push_back(new_pos);
        }

        if(!sg[level].sub_node[v].lambda.empty()){
            int lamda_list_len=sg[level].sub_node[v].lambda.size();
            //sg[level].sub_node[v].neighbour_s_node;
            for(int adj_sun=0;adj_sun<lamda_list_len;adj_sun++){
                int adj_sun_index=sg[level].sub_node[v].neighbour_s_node.at(adj_sun);
                lamda=sg[level].sub_node[v].lambda.at(adj_sun);
                adj_sun_pos=DPoint (sg[level].sub_node[adj_sun_index].m_x,sg[level].sub_node[adj_sun_index].m_y);
                new_pos= get_waggled_inbetween_position(sun_pos,adj_sun_pos,lamda);
                L.push_back(new_pos);
            }
        }
        DPoint barycenter= get_barycenter_position(L);
        sg[level].sub_node[v].m_x=barycenter.first;
        sg[level].sub_node[v].m_y=barycenter.second;
        sg[level].sub_node[v].placed=true;
        /*if(isnan(barycenter.first)){
            fprintf(stderr,"pm_nananannanananananan\n");
        }*/
        //fprintf(stderr,"pm_nodes    m_x=%f m_y=%f\n",sg[level].sub_node[v].m_x,sg[level].sub_node[v].m_y);
    }

}
void Multilevel::create_all_placement_sectors(int level,vector<subgraph> &sg){
    int num_of_level1=sg[level+1].sub_node.size();
   // cout<<"num of node"<<num_of_level1<<endl;
    for(int i=0;i<num_of_level1;i++){
        double angle_1=0,angle_2=0;
        vector<float> adj_pos_x,adj_pos_y;
        float v_high_mx=sg[level+1].sub_node[i].m_x;
        float v_high_my=sg[level+1].sub_node[i].m_y;
       // cout<<" i="<<i<<endl;
        for(int j=0;j<sg[level+1].edge[i].size();j++){
            int w_high=sg[level+1].edge[i][j].end;
            if(i==w_high) continue;
            adj_pos_x.push_back(sg[level+1].sub_node[w_high].m_x);
            adj_pos_y.push_back(sg[level+1].sub_node[w_high].m_y);
        }
       /* for (int j = 0; j < sg[level+1].sub_edge.size(); ++j) {
            int src,tgt;
            src=sg[level+1].sub_edge[j].source;
            tgt=sg[level+1].sub_edge[j].target;
            if(src==tgt) continue;

        }*/
        //cout<<" i="<<i<<endl;
        const float parallel_pos_x=v_high_mx+1;
        const float parallel_pos_y=v_high_my;
        if(adj_pos_x.size()==0){
            angle_2 = 2.0 * pi;
        } else if(adj_pos_x.size()==1){
            angle_1=angle(v_high_mx,v_high_my,parallel_pos_x,parallel_pos_y,adj_pos_x[0],adj_pos_y[0]);
           // fprintf(stderr,"angle=%f\n",angle_1);
            angle_2=angle_1+pi;
        }else{
            const int MAXX=10;
            int steps=1;
            int adj_len=adj_pos_x.size();
            for (int j = 0; j < adj_len; ++j) {
                double act_angle_1= angle(v_high_mx,v_high_my,parallel_pos_x,parallel_pos_y,adj_pos_x[j],adj_pos_y[j]);
               // fprintf(stderr,"act_angle_1=%.3f v_high_mx=%.3f v_high_my=%.3f parallel_pos_x=%.3f parallel_pos_y=%.3f adj_pos_x=%.3f adj_pos_y=%.3f\n",act_angle_1,v_high_mx,v_high_my,parallel_pos_x,parallel_pos_y,adj_pos_x[j],adj_pos_y[j]);
                double min_next_angle = std::numeric_limits<double>::max();
                for(int k=0;k<adj_len;k++){
                    if(k!=j){
                        double temp_angle= angle(v_high_mx,v_high_my,adj_pos_x[j],adj_pos_y[j],adj_pos_x[k],adj_pos_y[k]);
                        if(min_next_angle>temp_angle){
                            min_next_angle=temp_angle;
                        }
                    }
                }
                if(j==0||(min_next_angle>angle_2-angle_1)){
                    angle_1=act_angle_1;
                    angle_2=act_angle_1+min_next_angle;
                }
                steps++;
                if( steps>MAXX){
                    break;
                }

            }
            if(angle_1==angle_2){
                angle_2=angle_1+pi;
            }
            int sun_node=sg[level+1].sub_node[i].id;

            sg[level].sub_node[sun_node].angle_1=angle_1;
            sg[level].sub_node[sun_node].angle_2=angle_2;
        }

    }
    for(int i=0;i<sg[level].sub_node.size();i++){
        int ded_sun=sg[level].sub_node[i].dedicated_sun_node;
        sg[level].sub_node[i].angle_1=sg[level].sub_node[ded_sun].angle_1;
        sg[level].sub_node[i].angle_2=sg[level].sub_node[ded_sun].angle_2;
       //fprintf(stderr,"v=%d ded_sun=%d angle1=%f angle2=%f \n",i,ded_sun,sg[level].sub_node[i].angle_1,sg[level].sub_node[i].angle_2);
    }
}
float Multilevel::angle(float m_x,float m_y,float q_m_x,float q_m_y, float r_m_x,float r_m_y){
    const double dx1 = q_m_x - m_x, dy1 = q_m_y - m_y;
    const double dx2 = r_m_x - m_x, dy2 = r_m_y - m_y;

    // two vertices on the same place!
    if ((dx1 == 0 && dy1 == 0) || (dx2 == 0 && dy2 == 0)) {
        return 0.0;
    }

    double phi = std::atan2(dy2, dx2) - std::atan2(dy1, dx1);
    if (phi < 0) { phi += 2*pi; }

    return phi;
}

DPoint Multilevel::calculate_position(DPoint P,DPoint Q,float dist_P,float dist_Q){

        float dist_PQ= sqrt((P.first-Q.first)*(P.first-Q.first)+(P.second-Q.second)*(P.second-Q.second))+0.01;
        float lamba=(dist_P + (dist_PQ - dist_P - dist_Q)/2)/dist_PQ;
        /*if(isnan(dist_PQ)|| isnan(lamba)){
            fprintf(stderr,"dist_PQ or lamda error\n");
        }*/
        return get_waggled_inbetween_position(P,Q,lamba);

}
DPoint Multilevel::get_waggled_inbetween_position(DPoint s, DPoint t, double lambda)
{
    const double WAGGLEFACTOR = 0.05;
    const int BILLION = 1000000000;
    DPoint inbetween_point;
    inbetween_point.first = s.first + lambda*(t.first- s.first);
    inbetween_point.second = s.second + lambda*(t.second - s.second);
    float radius = WAGGLEFACTOR * norm(t,s);
   // srand(time(0));
    //int a=rand()%BILLION+1;
    //float rnd = float(a)/(float)(BILLION+2);//rand number in (0,1)
    double rnd = double(randomNumber(1,BILLION)+1)/(BILLION+2);//rand number in (0,1)
    float rand_radius =  radius * rnd;
    return create_random_pos(inbetween_point,rand_radius,0,2.0*pi);
}
float Multilevel::norm(DPoint P,DPoint Q){
    float norm= sqrt((P.first-Q.first)*(P.first-Q.first)+(P.second-Q.second)*(P.second-Q.second));
    return norm;
}
DPoint  Multilevel::create_random_pos(DPoint center,float radius,double angle_1,
                                            double angle_2)
{
    const int BILLION = 1000000000;
    DPoint  new_point;
    //srand(time(0));
   // int a=rand()%BILLION+1;
   // float rnd = float(a)/(float)(BILLION+2);//rand number in (0,1)
    double rnd = double(randomNumber(1,BILLION)+1)/(BILLION+2);//rand number in (0,1)
    double rnd_angle = angle_1 +(angle_2-angle_1)*rnd;
    double dx = cos(rnd_angle) * radius;
    double dy = sin(rnd_angle) * radius;
    new_point.first = center.first + dx ;
    new_point.second = center.second + dy;
    return new_point;
}
DPoint Multilevel::get_barycenter_position(vector<DPoint>& L){
    DPoint sum(0,0);
    DPoint barycenter;
    int len=L.size();
    for (int i = 0; i < len; ++i) {
        sum.first+=L[i].first;
        sum.second+=L[i].second;
    }
    barycenter.first=sum.first/len;
    barycenter.second=sum.second/len;
    return barycenter;
}
FMMMLayout::FMMMLayout()
{
    //fprintf(stderr,"FMMMLayout\n");
    initialize_all_options();
}
void FMMMLayout::call(graph &GA)
{
    vector<double> edgelength(GA.links.size(),1.0);
    call(GA,edgelength);
   // cout<<"111"<<endl;
}
void FMMMLayout::call(graph &GA, const vector<double> &edgeLength)
{
    //const Graph &G = GA.constGraph();
   // NodeArray<NodeAttributes> A(G);       //stores the attributes of the nodes (given by L)
    //EdgeArray<EdgeAttributes> E(G);       //stores the edge attributes of G
    graph G_reduced;                      //stores a undirected simple and loop-free copy of G
    //EdgeArray<EdgeAttributes> E_reduced;  //stores the edge attributes of G_reduced
    //NodeArray<NodeAttributes> A_reduced;  //stores the node attributes of G_reduced
    int numofnode=GA.nodes.size();
   // cout<<"1112222"<<endl;
    if(numofnode > 1)
    {
        //GA.clearAllBends();//all edges are straight-line
       // if(useHighLevelOptions())
            //update_low_level_options_due_to_high_level_options_settings();
        //import_NodeAttributes(G,GA,A);
        //import_EdgeAttributes(G,edgeLength,E);
        //init graph
        subgraph SG;
        for(int i=0;i<GA.nodes.size();i++){
            subNode sn;
            sn.id=i;
            sn.m_x=0;sn.m_y=0;
            sn.height=0;
            sn.weight=0;
            sn.v_higher_level_index=0;
            sn.v_lower_level_index=i;
            SG.sub_node.push_back(sn);
        }
        for(int i=0;i<GA.links.size();i++){
            adjLink al;
            al.source=GA.links[i].source;
            al.target=GA.links[i].target;
            if(edgeLength[i]>0) al.length=edgeLength[i];
            else al.length=1;

            SG.sub_edge.push_back(al);
        }
        max_integer_position = pow(2.0,40);
        //init_ind_ideal_edgelength(G,A,E);
        for(int i=0;i<SG.sub_node.size();i++){
            double w=SG.sub_node[i].weight/2.0;
            double h=SG.sub_node[i].height/2.0;
            radius.push_back(sqrt(w*w+ h*h));
        }
        for (int i = 0; i < SG.sub_edge.size(); ++i) {
            int source=SG.sub_edge[i].source;
            int target=SG.sub_edge[i].target;
            SG.sub_edge[i].length= SG.sub_edge[i].length*20.0+radius[source]+radius[target];
        }
        //cout<<"1111222111"<<endl;
        //make_simple_loopfree(G,A,E,G_reduced,A_reduced,E_reduced);
        //call_DIVIDE_ET_IMPERA_step(G_reduced,A_reduced,E_reduced);
        call_MULTILEVEL_step_for_subGraph(SG);
       // pack_subGraph_drawings (A,G_sub,A_sub);
        //delete_all_subGraphs(G_sub,A_sub,E_sub);
        //adjust_positions(G_reduced, A_reduced);


        //export_NodeAttributes(G_reduced,A_reduced,GA);
        vector<int> flag(GA.nodes.size(),-1);
        for (int i = 0; i < SG.sub_node.size(); ++i) {
            GA.nodes[SG.sub_node[i].id].x=SG.sub_node[i].m_x;
            GA.nodes[SG.sub_node[i].id].y=SG.sub_node[i].m_y;
            flag[SG.sub_node[i].id]=1;
          // GA.nodes[i].x=SG.sub_node[i].m_x;
            //GA.nodes[i].y=SG.sub_node[i].m_y;
        }
        for (int i = 0; i < flag.size(); ++i) {
            if(flag[i]==-1){
                fprintf(stderr,"the node is not placed\n");
                //cout<<"the node all are placed"<<endl;
            }
        }
        fprintf(stderr,"the node all are placed\n");
        //cout<<"the node all are placed"<<endl;
    }
    else //trivial cases
    {
       GA.nodes[0].x=0;
       GA.nodes[0].y=0;
    }
}
void FMMMLayout::call_MULTILEVEL_step_for_subGraph(
        subgraph& SG){
    int max_level = 30;
    if (m_singleLevel) m_minGraphSize = SG.sub_node.size();
    vector<subgraph> sg(max_level+1);
    SG.numofnode=SG.sub_node.size();
    SG.numofedge=SG.sub_edge.size();
    mul.create_multilevel_representations(SG,sg,max_level,randSeed(),
                                           minGraphSize(),randomTries());
    for(int i = max_level;i >= 0;i--)
    {
        cout<<"level"<<i<<endl;
        if(i == max_level){
            create_initial_placement(sg,i);
        }else
        {

            mul.find_initial_placement_for_level(i,sg);
            update_boxlength_and_cornercoordinate(sg,i);
        }
       // cout<<"level"<<i<<endl;
        call_FORCE_CALCULATION_step(sg, i,max_level);
    }
    for (int i = 0; i < SG.sub_node.size(); ++i) {
        SG.sub_node[i].m_x=sg[0].sub_node[i].m_x;
        SG.sub_node[i].m_y=sg[0].sub_node[i].m_y;
    }
   // Mult.delete_multilevel_representations(G_mult_ptr,A_mult_ptr,E_mult_ptr,max_level);
}
void FMMMLayout::create_initial_placement (vector<subgraph>& sg,
                               int act_level){
    init_boxlength_and_cornercoordinate(sg,act_level);
    mul.setSeed((unsigned int) time(nullptr));
    const int BILLION = 1000000000;

    for (int v=0;v<sg[act_level].sub_node.size();v++) {
        DPoint rndp;
        rndp.first = double(mul.randomNumber(0, BILLION)) / BILLION; //rand_x in [0,1]
        rndp.second = double(mul.randomNumber(0, BILLION)) / BILLION; //rand_y in [0,1]
        sg[act_level].sub_node[v].m_x=(rndp.first*(boxlength - 2) + 1);
        sg[act_level].sub_node[v].m_y=(rndp.second*(boxlength - 2) + 1);
        //fprintf(stderr,"v=%d m_x=%f m_y=%f\n",v,sg[act_level].sub_node[v].m_x,sg[act_level].sub_node[v].m_y);
    }

    update_boxlength_and_cornercoordinate(sg,act_level);
}
void FMMMLayout::init_boxlength_and_cornercoordinate(vector<subgraph> &sg,int level){
    const double MIN_NODE_SIZE = 10;
    const double BOX_SCALING_FACTOR = 1.1;

    double w = 0, h = 0;
    for (int v ;v<sg[level].sub_node.size();v++)
    {
        w += max(sg[level].sub_node[v].weight, MIN_NODE_SIZE);
        h += max(sg[level].sub_node[v].height, MIN_NODE_SIZE);
    }

    boxlength = ceil(max(w, h) * BOX_SCALING_FACTOR);
    //fprintf(stderr,"boxlen=%f\n",boxlength);
//down left corner of comp. box is the origin
    down_left_corner.first = 0;
    down_left_corner.second = 0;
}

void FMMMLayout::update_boxlength_and_cornercoordinate(vector<subgraph> &sg,int level){
    int node_of_level=sg[level].sub_node.size();

    float mid_x=sg[level].sub_node[0].m_x;
    float mid_y=sg[level].sub_node[0].m_y;


    double xmin, xmax, ymin, ymax;
    xmin = xmax = mid_x;
    ymin = ymax = mid_y;

    for (int i = 1; i < node_of_level; ++i) {
        mid_x=sg[level].sub_node[i].m_x;
        mid_y=sg[level].sub_node[i].m_y;
        if (mid_x < xmin)
            xmin = mid_x;
        if (mid_x > xmax)
            xmax = mid_x;
        if (mid_y < ymin)
            ymin = mid_y;
        if (mid_y > ymax)
            ymax = mid_y;
    }

    //set down_left_corner and boxlength

    down_left_corner.first = floor(xmin - 1);
    down_left_corner.second = floor(ymin - 1);
    boxlength = ceil(max(ymax - ymin, xmax - xmin) *1.01 + 2);

    //exception handling: all nodes have same x and y coordinate
    if (boxlength <= 2)
    {
        boxlength = node_of_level* 20;
        down_left_corner.first = floor(xmin) - (boxlength / 2);
        down_left_corner.second = floor(ymin) - (boxlength / 2);
    }
    //fprintf(stderr,"boxlen=%f down_left_corner.mx=%f down_left_corner.my=%f\n",boxlength,down_left_corner.first,down_left_corner.second);
    //export the boxlength and down_left_corner values to the rep. calc. classes
    FR.update_boxlength_and_cornercoordinate(boxlength, down_left_corner);
    //fprintf(stderr,"boxlen=%f down_left_corner.mx=%f down_left_corner.my=%f\n",boxlength,down_left_corner.first,down_left_corner.second);

}
void FMMMLayout::initialize_all_options()
{
    fprintf(stderr,"FMMMLayout\n");
    //setting high level options
    useHighLevelOptions(false);
  //  pageFormat(FMMMOptions::PageFormatType::Square);
  //  unitEdgeLength(LayoutStandards::defaultNodeSeparation());
    newInitialPlacement(false);
   // qualityVersusSpeed(FMMMOptions::QualityVsSpeed::BeautifulAndFast);

    //setting low level options
    //setting general options
    randSeed(100);
   // fprintf(stderr,"randseed=%d\n",randSeed());
    //edgeLengthMeasurement(FMMMOptions::EdgeLengthMeasurement::BoundingCircle);
   // allowedPositions(FMMMOptions::AllowedPositions::Integer);
    maxIntPosExponent(40);
   // fprintf(stderr,"maxIntPosExponent=%d\n",maxIntPosExponent());
    //setting options for the divide et impera step
    pageRatio(1.0);
    stepsForRotatingComponents(10);
    //tipOverCCs(FMMMOptions::TipOver::NoGrowingRow);
   // minDistCC(LayoutStandards::defaultCCSeparation());
   // presortCCs(FMMMOptions::PreSort::DecreasingHeight);

    //setting options for the multilevel step
    minGraphSize(50);
    //galaxyChoice(FMMMOptions::GalaxyChoice::NonUniformProbLowerMass);
    randomTries(20);
   // maxIterChange(FMMMOptions::MaxIterChange::LinearlyDecreasing);
    maxIterFactor(10);
    //initialPlacementMult(FMMMOptions::InitialPlacementMult::Advanced);
    m_singleLevel = false;

    //setting options for the force calculation step
   // forceModel(FMMMOptions::ForceModel::New);
    springStrength(1);
    repForcesStrength(1);
    //repulsiveForcesCalculation(FMMMOptions::RepulsiveForcesMethod::NMM);
    //stopCriterion(FMMMOptions::StopCriterion::FixedIterationsOrThreshold);
    threshold(0.01);
    fixedIterations(30);
    forceScalingFactor(0.05);
    coolTemperature(false);
    coolValue(0.99);
    //initialPlacementForces(FMMMOptions::InitialPlacementForces::RandomRandIterNr);

    //setting options for postprocessing
    resizeDrawing(true);
    resizingScalar(1);
    fineTuningIterations(40);
    fineTuneScalar(0.2);
    adjustPostRepStrengthDynamically(true);
    postSpringStrength(2.0);
    postStrengthOfRepForces(0.01);

    //setting options for different repulsive force calculation methods
    frGridQuotient(2);
   // nmTreeConstruction(FMMMOptions::ReducedTreeConstruction::SubtreeBySubtree);
    //nmSmallCell(FMMMOptions::SmallestCellFinding::Iteratively);
    nmParticlesInLeaves(25);
    nmPrecision(4);
}
void FMMMLayout::call_FORCE_CALCULATION_step(vector<subgraph> &sg,int act_level,int max_level){
    int num_of_node=sg[act_level].sub_node.size();
    cout<<"force calculation"<<endl;
    if(num_of_node>1){
        int iter=1;
        int max_mult_iter=mul.get_max_mult_iter(act_level,max_level,num_of_node);
        //cout<<"max_iter"<<max_mult_iter<<endl;
        double threshold=0.01;
        double actforcevectorlength = threshold + 1;

        vector<DPoint> F_rep(num_of_node,DPoint(0,0));
        vector<DPoint> F_attr(num_of_node,DPoint(0,0));
        vector<DPoint> F(num_of_node,DPoint(0,0));
        vector<DPoint> last_node_movement(num_of_node,DPoint(0,0));
        set_average_ideal_edgelength(sg,act_level);
      //  fprintf(stderr,"average_ideal_edgelength=%f actforcevectorlength=%f \n",average_ideal_edgelength,actforcevectorlength);
        make_initialisations_for_rep_calc_classes(sg,act_level);
        //cout<<"force calculation"<<endl;
        int running_time=0;
        while (running(iter, max_mult_iter, actforcevectorlength)) {
           // fprintf(stderr,"running_time=%d\n",running_time);
            running_time++;
            calculate_forces(sg,act_level,F,F_attr,F_rep,last_node_movement,iter,0);
            actforcevectorlength = get_average_forcevector_length(sg,act_level,F);
            iter++;
        }

        if(act_level == 0)
            call_POSTPROCESSING_step(sg,0,F,F_attr,F_rep,last_node_movement);

       // deallocate_memory_for_rep_calc_classes();
    }
}
bool FMMMLayout::running(int iter, int max_mult_iter, double actforcevectorlength)
{
    const int ITERBOUND = 10000;
    return iter <= max_mult_iter && actforcevectorlength >= threshold();
    return false;
}
void FMMMLayout::call_POSTPROCESSING_step(
        vector<subgraph> &sg,
        int act_level,
        vector<DPoint>& F,
        vector<DPoint>& F_attr,
        vector<DPoint>& F_rep,
        vector<DPoint>& last_node_movement){
    for(int i = 1; i<= 10; i++)
        calculate_forces(sg,act_level,F,F_attr,F_rep,last_node_movement,i,1);

    if(resizeDrawing())
    {
      //  fprintf(stderr,"resizedrawing111111\n");
        adapt_drawing_to_ideal_average_edgelength(sg,act_level);
        update_boxlength_and_cornercoordinate(sg,act_level);
    }

    for(int i = 1; i<= fineTuningIterations(); i++)
        calculate_forces(sg,act_level,F,F_attr,F_rep,last_node_movement,i,2);

    if(resizeDrawing()){
       // fprintf(stderr,"resizedrawing2222\n");
        adapt_drawing_to_ideal_average_edgelength(sg,act_level);
    }

}
void FMMMLayout::adapt_drawing_to_ideal_average_edgelength(vector<subgraph>& sg,int act_level){
    double sum_real_edgelength = 0;
    double sum_ideal_edgelength = 0;
    double real_edgelength = 0;
    int num_of_edge=sg[act_level].sub_edge.size();
    for (int i = 0; i < num_of_edge; ++i)
    {
        sum_ideal_edgelength += sg[act_level].sub_edge[i].length;
        int src=sg[act_level].sub_edge[i].source;
        int tgt=sg[act_level].sub_edge[i].target;
        float mvx=sg[act_level].sub_node[src].m_x-sg[act_level].sub_node[tgt].m_x;
        float mvy=sg[act_level].sub_node[src].m_y-sg[act_level].sub_node[tgt].m_y;
        sum_real_edgelength += sqrt(mvx*mvx+mvy*mvy);
    }

    double area_scaling_factor;
    if (sum_real_edgelength == 0) // very very unlikely case
        area_scaling_factor = 1;
    else
        area_scaling_factor = sum_ideal_edgelength / sum_real_edgelength;
   // fprintf(stderr,"area_scaling_factor=%f \n",area_scaling_factor);
    DPoint new_pos;
    int num_of_node=sg[act_level].sub_node.size();
    float resizingScalar=1.0;
    for (int i = 0; i < num_of_node; ++i) {
        new_pos.first = resizingScalar * area_scaling_factor * sg[act_level].sub_node[i].m_x;
        new_pos.second = resizingScalar * area_scaling_factor * sg[act_level].sub_node[i].m_y;
        sg[act_level].sub_node[i].m_x=new_pos.first;
        sg[act_level].sub_node[i].m_y=new_pos.second;
    }



}
double FMMMLayout::get_average_forcevector_length (vector<subgraph>& sg,int act_level, vector<DPoint>& F){
    double lengthsum = 0;
    int num_of_node=sg[act_level].sub_node.size();
    for (int i = 0; i < num_of_node; ++i) {
        lengthsum += sqrt(F[i].second*F[i].second+F[i].first*F[i].first);
    }

    lengthsum /= num_of_node;
    return lengthsum;
}
void FMMMLayout::calculate_forces(
        vector<subgraph> &sg,int act_level,
        vector<DPoint>& F,
        vector<DPoint>& F_attr,
        vector<DPoint>& F_rep,
        vector<DPoint>& last_node_movement,
        int iter,
        int fine_tuning_step)
{
    //adjust_positions( sg,act_level);
    calculate_attractive_forces(sg,act_level,F_attr);
    /*for (int i = 0; i <sg[act_level].sub_node.size(); ++i) {
        cout<<sg[act_level].sub_node[i].m_x<< " "<<sg[act_level].sub_node[i].m_y<<",";
    }*/
    //cout<<"before move"<<sg[act_level].sub_node[0].m_x<< " "<<sg[act_level].sub_node[0].m_y<<endl;
    //cout<<"Fa"<<F_attr[0].first<< "  "<<F_attr[0].second<<endl;
    calculate_repulsive_forces(sg,act_level,F_rep);
    //cout<<"Fr"<<F_rep[0].first<< " "<<F_rep[0].second<<endl;
    add_attr_rep_forces(sg,act_level,F_attr,F_rep,F,iter,fine_tuning_step);
   //cout<<"add"<<F[0].first<< " "<<F[0].second<<endl;
    prevent_oscillations(sg,act_level,F,last_node_movement,iter);
   // cout<<"prevent_oscillations"<<F[0].first<< " "<<F[0].second<<endl;
    move_nodes(sg,act_level,F);
    //cout<<"move"<<sg[act_level].sub_node[0].m_x<< " "<<sg[act_level].sub_node[0].m_y<<endl;
    update_boxlength_and_cornercoordinate(sg,act_level);
}
void FMMMLayout::adjust_positions( vector<subgraph>& sg,int act_level){


}
void FMMMLayout:: calculate_attractive_forces(vector<subgraph>& sg,int act_level,vector<DPoint>& F_attr){
    DPoint f_u;
    DPoint nullpoint (0,0);

    //initialisation
   // init_F(G,F_attr);
    int num_of_node=sg[act_level].sub_node.size();
    int num_of_edge=sg[act_level].sub_edge.size();
    for (int i = 0; i < num_of_node; ++i) {
        F_attr[i] = nullpoint;
    }
    //calculation
    for (int i = 0; i < num_of_edge; ++i) {
        int u = sg[act_level].sub_edge[i].source;
        int v =  sg[act_level].sub_edge[i].target;

        DPoint vector_v_minus_u ;
        vector_v_minus_u.first= sg[act_level].sub_node[v].m_x-sg[act_level].sub_node[u].m_x;
        vector_v_minus_u.second= sg[act_level].sub_node[v].m_y-sg[act_level].sub_node[u].m_y;
       // if(act_level==0)
        //    fprintf(stderr," u_mx=%f u_my=%f v_mx=%f v_my=%f vector_v_minus_u.first=%f vector_v_minus_u.second=%f  \n",sg[act_level].sub_node[u].m_x,sg[act_level].sub_node[u].m_y,sg[act_level].sub_node[v].m_x,sg[act_level].sub_node[v].m_y,vector_v_minus_u.first,vector_v_minus_u.second);
        double norm_v_minus_u = sqrt(vector_v_minus_u.first*vector_v_minus_u.first+ vector_v_minus_u.second*vector_v_minus_u.second);
        if(vector_v_minus_u == nullpoint){
            f_u = nullpoint;
            //fprintf(stderr,"-----------------f_near_machine_precision src=%d tgt=%d\n",u,v);
        }

        else if(!mul.f_near_machine_precision(norm_v_minus_u,f_u))
        {
           // fprintf(stderr,"!f_near_machine_precision src=%d tgt=%d\n",u,v);
            double scalar = f_attr_scalar(norm_v_minus_u,sg[act_level].sub_edge[i].length)/norm_v_minus_u;
           //if(act_level==0)
           //  fprintf(stderr," src=%d tgt=%d scalar=%f norm_v_minus_u=%f length=%f \n",u,v,scalar,norm_v_minus_u,sg[act_level].sub_edge[i].length);
            f_u.first = scalar * vector_v_minus_u.first;
            f_u.second = scalar * vector_v_minus_u.second;
        }

        F_attr[v].first = F_attr[v].first - f_u.first;
        F_attr[v].second = F_attr[v].second - f_u.second;
        F_attr[u].first = F_attr[u].first + f_u.first;
        F_attr[u].second = F_attr[u].second + f_u.second;
    }

}
double FMMMLayout::f_attr_scalar(double d, double ind_ideal_edge_length)
{
    double s(0);
    const double c =  std::log2(d/ind_ideal_edge_length);
    if (d > 0)
        s =  c * d * d /(ind_ideal_edge_length * ind_ideal_edge_length * ind_ideal_edge_length);
    else
        s = -1e10;

    return s;
}

double Multilevel::random_precision_number(double shift)
{
    const int BILLION = 1000000000;
    double rand = shift + double(randomNumber(1, BILLION) + 1) / (BILLION + 2);
    return randomNumber(0, 1) == 0 ? rand : -rand;
}
bool Multilevel::f_near_machine_precision(double distance, DPoint& force)
{
    const double POS_BIG_LIMIT = POS_BIG_DOUBLE * 1e-190;
    const double POS_SMALL_LIMIT = POS_SMALL_DOUBLE * 1e190;

    if (distance < POS_SMALL_LIMIT) {
        force = DPoint(POS_SMALL_LIMIT * random_precision_number(1),
                       POS_SMALL_LIMIT * random_precision_number(1));
        return true;
    } else if (distance > POS_BIG_LIMIT) {
        force = DPoint(POS_BIG_LIMIT * random_precision_number(0),
                       POS_BIG_LIMIT * random_precision_number(0));
        return true;
    }
    return false;
}
double Multilevel::randomDouble(double low, double high)
{


    std::uniform_real_distribution<> dist(low,high);

#ifndef OGDF_MEMORY_POOL_NTS
    std::lock_guard<std::mutex> guard(s_randomMutex);
#endif
    return dist(s_random);
}

void FMMMLayout::calculate_repulsive_forces(vector<subgraph>& sg,int act_level,vector<DPoint>& F_rep){
    FR.calculate_exact_repulsive_forces(sg,act_level,F_rep);
}
void  FMMMLayout::add_attr_rep_forces(vector<subgraph>& sg,int act_level,vector<DPoint>& F_attr,vector<DPoint>& F_rep,vector<DPoint>& F,int iter,int fine_tuning_step){
    DPoint nullpoint(0, 0);
    int num_of_node=sg[act_level].sub_node.size();
    //set cool_factor
    if (!coolTemperature())
        cool_factor = 1.0;
    else if (coolTemperature() && fine_tuning_step == 0)
    {
        if (iter == 1)
            cool_factor = coolValue();
        else
            cool_factor *= coolValue();
    }

    if (fine_tuning_step == 1)
        cool_factor /= 10.0; //decrease the temperature rapidly
    else if (fine_tuning_step == 2)
    {
        if (iter <= fineTuningIterations() - 5)
            cool_factor = fineTuneScalar(); //decrease the temperature rapidly
        else
            cool_factor = (fineTuneScalar() / 10.0);
    }

    //set the values for the spring strength and strength of the rep. force field
    double act_spring_strength, act_rep_force_strength;
    if (fine_tuning_step <= 1)//usual case
    {
        act_spring_strength = springStrength();
        act_rep_force_strength = repForcesStrength();
    }
    else if (!adjustPostRepStrengthDynamically())
    {
        act_spring_strength = postSpringStrength();
        act_rep_force_strength = postStrengthOfRepForces();
    }
    else //adjustPostRepStrengthDynamically())
    {
        act_spring_strength = postSpringStrength();
        act_rep_force_strength = get_post_rep_force_strength(num_of_node);
    }
    //fprintf(stderr,"cool_factor=%f coolTemperature=%d fine_tuning_step=%d coolValue=%f fineTuningIterations()=%d fineTuneScalar=%f springStrength=%f repForcesStrength=%f adjustPostRepStrengthDynamically=%d\n",
     //       cool_factor,coolTemperature(),fine_tuning_step,coolValue(),fineTuningIterations(),fineTuneScalar(),springStrength(),repForcesStrength(),adjustPostRepStrengthDynamically());
    // fprintf(stderr,"act_spring_strength=%f act_rep_force_strength=%f \n",act_spring_strength,act_rep_force_strength);
    for (int v = 0; v < num_of_node; ++v) {

        DPoint f;
        f.first = act_spring_strength * F_attr[v].first + act_rep_force_strength * F_rep[v].first;
        f.second = act_spring_strength * F_attr[v].second + act_rep_force_strength * F_rep[v].second;
        f.first = average_ideal_edgelength * average_ideal_edgelength * f.first;
        f.second = average_ideal_edgelength * average_ideal_edgelength * f.second;

        double norm_f = sqrt(f.first*f.first+f.second*f.second);

        DPoint force;
        if (f == nullpoint)
            force = nullpoint;
        else if (mul.f_near_machine_precision(norm_f, force))
            restrict_force_to_comp_box(force);
        else
        {
            double scalar = min(norm_f * cool_factor * forceScalingFactor(), max_radius(iter)) / norm_f;
            force.first = scalar * f.first;
            force.second = scalar * f.second;
        }
        F[v] = force;
    }
}
float Multilevel:: angle( DPoint q, DPoint r,DPoint a) {
    const double dx1 = q.first - a.first, dy1 = q.second - a.second;
    const double dx2 = r.first - a.first, dy2 = r.second - a.second;

    // two vertices on the same place!
    if ((dx1 == 0 && dy1 == 0) || (dx2 == 0 && dy2 == 0)) {
        return 0.0;
    }

    double phi = std::atan2(dy2, dx2) - std::atan2(dy1, dx1);
    if (phi < 0) { phi += 2*pi; }

    return phi;
}
void  FMMMLayout::prevent_oscillations(vector<subgraph>& sg,int act_level,vector<DPoint>& F,vector<DPoint>&last_node_movement,int iter){
    const double pi_times_1_over_6 = 0.52359878;
    const double factors[] = {
            2.0, 2.0, 1.5, 1.0, 0.66666666, 0.5, 0.33333333,
            0.33333333, 0.5, 0.66666666, 1.0, 1.5, 2.0, 2.0
    };
    const DPoint nullpoint(0, 0);

    if (iter > 1) { // usual case
        int num_of_node=sg[act_level].sub_node.size();
        for (int v = 0; v < num_of_node; ++v) {
            const DPoint force_new(F[v]);
            const DPoint force_old(last_node_movement[v]);
            const double norm_new = sqrt(F[v].first*F[v].first+F[v].second*F[v].second);
            const double norm_old = sqrt(last_node_movement[v].first*last_node_movement[v].first+last_node_movement[v].second*last_node_movement[v].second);
            if (norm_new > 0 && norm_old > 0) {
                const double fi = mul.angle(force_old, force_new,nullpoint);
                const double factor = factors[int(std::ceil(fi / pi_times_1_over_6))];
                const double quot = norm_old * factor / norm_new;
                if (quot < 1.0) {
                    F[v].first = quot * F[v].first;
                    F[v].second = quot * F[v].second;
                }
            }
            last_node_movement[v] = F[v];
        }
    }
    else if (iter == 1)
        init_last_node_movement(sg,act_level,F,last_node_movement);
}
void FMMMLayout::init_last_node_movement(
        vector<subgraph>& sg,
        int act_level,
        vector<DPoint>& F,
        vector<DPoint>& last_node_movement){
    int num_of_node=sg[act_level].sub_node.size();
    for (int v = 0; v < num_of_node; ++v) {
        last_node_movement[v]= F[v];
    }

}
double angle(DPoint a,DPoint q, DPoint r) {
    const double dx1 = q.first - a.first, dy1 = q.second - a.second;
    const double dx2 = r.first - a.first, dy2 = r.second - a.second;

    // two vertices on the same place!
    if ((dx1 == 0 && dy1 == 0) || (dx2 == 0 && dy2 == 0)) {
        return 0.0;
    }

    double phi = std::atan2(dy2, dx2) - std::atan2(dy1, dx1);
    if (phi < 0) { phi += 2*pi; }

    return phi;
}
void FMMMLayout::move_nodes(vector<subgraph>& sg,int act_level,vector<DPoint>& F){
    int num_of_node=sg[act_level].sub_node.size();
    for (int v = 0; v < num_of_node; ++v) {
        sg[act_level].sub_node[v].m_x +=F[v].first;
        sg[act_level].sub_node[v].m_y +=F[v].second;
    }

}

int Multilevel::get_max_mult_iter(int act_level, int max_level, int node_nr){
    int maxIterFactor=10,fixedIterations=60;
    int iter;
    if(max_level==0){
        iter = maxIterFactor* fixedIterations;
    }else{
        iter = fixedIterations +
               int((double(act_level)/double(max_level)) *
                   (maxIterFactor - 1) * fixedIterations);
    }
    return iter;
    }

void FMMMLayout::set_average_ideal_edgelength(vector<subgraph> &sg,int act_level){
    int num_of_edge=sg[act_level].sub_edge.size();
    if(num_of_edge>0){
        double averagelength = 0;
        for (int i = 0; i < num_of_edge; ++i) {
            averagelength+=sg[act_level].sub_edge[i].length;
        }
        average_ideal_edgelength=averagelength/num_of_edge;
    }else{
        average_ideal_edgelength=50;
    }
}
void FMMMLayout::make_initialisations_for_rep_calc_classes(vector<subgraph>& sg,
        int act_level){
    int MIN_NODE_NUMBER=175;
   //fprintf(stderr,"boxlength=%f down_left_corner_mx=%f down_left_corner_my=%f frGridQuotient=%f\n",boxlength,down_left_corner.first,down_left_corner.second,frGridQuotient());
    //float frGridQuotient=1;
    frGridQuotient(2);
    FR.make_initialisations(boxlength, down_left_corner, 2);
    //fprintf(stderr,"boxlength=%f down_left_corner_mx=%f down_left_corner_my=%f frGridQuotient=%f\n",boxlength,down_left_corner.first,down_left_corner.second,frGridQuotient());
    /*if (sg[act_level].sub_node.size() >= MIN_NODE_NUMBER) { // using_NMM
        //using_NMM = true; //indicate that NMM is used for force calculation

        particles_in_leaves(p_i_l);
        precision(p);
        tree_construction_way(t_c_w);
        find_sm_cell(f_s_c);
        down_left_corner = d_l_c; //Export this two values from FMMM
        boxlength = bl;
        init_binko(2* precision());
    } else { // use exact method
        using_NMM = false; //indicate that exact method is used for force calculation
        ExactMethod.make_initialisations(bl,d_l_c,0);
    }*/

}