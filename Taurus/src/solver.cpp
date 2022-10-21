#include "solver.h"

int str2int(string s) {
    stringstream stream(s);
    int i;
    stream >> i;
    stream.clear();
    return i;
}
float str2float(string s) {
    stringstream stream(s);
    float i;
    stream >> i;
    stream.clear();
    return i;
}


void QuadForce(Quad* qt, float& fx, float& fy, float posx, float posy, float theta, float alpha, int nodeindex,float dp,int layer=0);
vector<int> split(const string& str, const string& delim) {
    vector<int> res;
    if ("" == str) return res;
    char * strs = new char[str.length() + 1];
    strcpy(strs, str.c_str());

    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while (p) {
        string s = p;
        int y = atoi(s.c_str());
        res.push_back(y);
        p = strtok(NULL, d);
    }
    return res;
}
long
rk_long(rk_state *state)
{
    return rk_ulong(state) >> 1;
}

unsigned long
rk_ulong(rk_state *state)
{
#if ULONG_MAX <= 0xffffffffUL
    return rk_random(state);
#else
    return (rk_random(state) << 32) | (rk_random(state));
#endif
}
unsigned long
rk_random(rk_state *state)
{

    unsigned long y;

    if (state->pos == RK_STATE_LEN) {
        int i;

        for (i = 0; i < 624 - 397; i++) {
            y = (state->key[i] & 0x80000000UL) | (state->key[i+1] & 0x7fffffffUL);
            state->key[i] = state->key[i+397] ^ (y>>1) ^ (-(y & 1) & 0x9908b0dfUL);
        }
        for (; i < 624 - 1; i++) {
            y = (state->key[i] & 0x80000000UL) | (state->key[i+1] & 0x7fffffffUL);
            state->key[i] = state->key[i+(397-624)] ^ (y>>1) ^ (-(y & 1) & 0x9908b0dfUL);
        }
        y = (state->key[624 - 1] & 0x80000000UL) | (state->key[0] & 0x7fffffffUL);
        state->key[624 - 1] = state->key[397 - 1] ^ (y >> 1) ^ (-(y & 1) & 0x9908b0dfUL);

        state->pos = 0;
    }
    y = state->key[state->pos++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
unsigned long
rk_interval(unsigned long max, rk_state *state)
{
    unsigned long mask = max, value;

    if (max == 0) {
        return 0;
    }
    /* Smallest bit mask >= max */
    mask |= mask >> 1;
    mask |= mask >> 2;
    mask |= mask >> 4;
    mask |= mask >> 8;
    mask |= mask >> 16;
#if ULONG_MAX > 0xffffffffUL
    mask |= mask >> 32;
#endif

    /* Search a random value in [0..mask] <= max */
#if ULONG_MAX > 0xffffffffUL
    if (max <= 0xffffffffUL) {
        while ((value = (rk_random(state) & mask)) > max);
    }
    else {
        while ((value = (rk_ulong(state) & mask)) > max);
    }
#else
    while ((value = (rk_ulong(state) & mask)) > max);
#endif
    return value;
}
void
rk_seed(unsigned long seed, rk_state *state)
{
    int pos;
    seed &= 0xffffffffUL;

    /* Knuth's PRNG as used in the Mersenne Twister reference implementation */
    for (pos = 0; pos < RK_STATE_LEN; pos++) {
        state->key[pos] = seed;
        seed = (1812433253UL * (seed ^ (seed >> 30)) + pos + 1) & 0xffffffffUL;
    }
    state->pos = RK_STATE_LEN;
    state->gauss = 0;
    state->has_gauss = 0;
    state->has_binomial = 0;
}
graph::graph() {
    n=0;m=0;
}

/*void graph::readGraph(string filename){
    ifstream input(filename);
    string s;
    string b="break";
    while (getline(input, s)) {
        if (s.at(0) != '%') {
            break;
        }
    }
    cout<<s<<endl;
    vector<int> di = split(s, " ");
    if(di.size()==3){
        n = di[1];
        m = di[2];

        //vector<Link> links = vector<Link>(m);
        for(int i=0;i<m;i++){
            if(!getline(input,s)) cout<<"file read error!"<<endl;
            vector<int> di = split(s, " ");
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
            vector<int> di1 = split(s, " ");
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
    cout<<n<<" "<<m<<endl;
}*/

void graph::initgraph(string filename) {
    ifstream input(filename);
    string s;
    string b="break";
    while (getline(input, s)) {
        if (s.at(0) != '%') {
            break;
        }
    }
    vector<int> I,J;

    vector<int> di = split(s, " ");
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
    vector<std::set<int> > undirected(n);

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
    cout<<"build undirected graph n="<<n<<" m="<<m<<endl;
}

void BFS(int V,int E,vector<Edge> e[MAX],int st,vector<int> &ind,vector<int>& count_com){
    //int visited[V];
    int *visited=new int[V];
    memset(visited,0,sizeof(visited));
    queue<int> q1;
    q1.push(st);
    //count_com[ind[st]]++;
    while(!q1.empty()){
        int u=q1.front();
        q1.pop();
        if(visited[u]) continue;
        visited[u]=1;
        ind[u]=ind[st];
        //cout<<"indu="<<ind[u]<<"  indst"<<ind[st]<<endl;
        count_com[ind[u]]++;
        for(int i=0;i<e[u].size();i++){
            int v=e[u][i].end;
            if(!visited[v]){
                q1.push(v);

            }
        }
    }
}
void graph::readGraphMtx(string filename){
    ifstream input(filename);
    string s;
    string b="break";
    int src,tgt,num_n=0,num_m=0;
    while(getline(input, s)){
        if(s==b){break;}
        num_m++;
        vector<int> di = split(s, " ");
        Link temp;
        src=di[0];
        tgt=di[1];
        temp.source=src;
        temp.target=tgt;
        links.push_back(temp);
        num_n=max(num_n,src);
        num_n=max(num_n,tgt);
    }

    n=num_n+1;
    m=num_m;

    vector<int> ind(n,-1);
    vector<int> count_com(n,0);

    //max component
    //vector<Edge> e[n];
    vector<Edge> *e=new vector<Edge>[n];
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
            //cout<<"indi=="<<ind[i]<<endl;
            BFS(n,m,e,i,ind,count_com);
            count++;
        }
    }
    if(count>1){
        //cout<<"The graph is disconnected! I am taking the largest component!"<<endl;

        int max_com=0,max=0;
        for(int i=0;i<count;i++){
            int temp=count_com[i];
            if(temp>max){max=temp;max_com=i;}
        }
        //cout<<"The component has="<<count_com[max_com]<<" vertex!"<<endl;

        count=0;
        map<int,int> map_ind;
        for(int i=0;i<n;i++){
            if(ind[i]!=max_com){
                count++;
                map_ind.insert(pair<int,int>(i,-1));
            }else{
                map_ind.insert(pair<int,int>(i,i-count));
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

    if(nodes.size()!=n){nodes.resize(n);}
    if(p>n) p=n;
    int I[m];
    int J[m];
    for (int i = 0; i < m; i++) {
        I[i]=links[i].source;
        J[i]=links[i].target;
    }
    vector<vector<int> > g1 = build_graph_unweighted(n, m, I, J);
    vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 48);
    set<int> pset;
    for(int i=0;i<closest_pivots.size();i++){
        pset.insert(closest_pivots[i]);
    }
    vector<int> pivot;

    for(set<int>::iterator it=pset.begin();it!=pset.end();it++){
        pivot.push_back(*it);
    }

    cout<<"compute pivotMDS"<<endl;
    PivotMDS(nodes,pivot.size(),n,sp,pivot);
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
    }
}
void graph::initNode(int n) {
    nodes.resize(n);
}
void graph::initEdge(int m) {
    links.resize(m);
}
void Dij(int V,int E,vector<Edge> e[MAX],int st,vector<float> &dist){
    priority_queue<Pair,vector<Pair>,greater<Pair>> q;
    int visited[V];
    memset(visited,0,sizeof(visited));
    for(int i=0;i<V;i++) dist[i]=MAX;
    dist[st]=0;
    q.push(make_pair(0,st));
    while(!q.empty()){
        Pair t=q.top();
        q.pop();
        if(visited[t.second]) continue;
        visited[t.second]=1;

        for(int i=0;i<e[t.second].size();i++){
            int son=e[t.second][i].end;
            if(dist[son]>dist[t.second]+e[t.second][i].w){
                dist[son]=dist[t.second]+e[t.second][i].w;
                if(!visited[son])
                    q.push(make_pair(dist[son],son));
            }
        }
    }
}
void maxmin_bfs_unweighted(const vector<vector<int> >& graph, const int p, vector<int>& mins, vector<int>& argmins)
{
    int n = graph.size();
    std::queue<int> q;
    vector<int> d(n, -1);

    q.push(p);
    d[p] = 0;
    while (!q.empty())
    {
        int current = q.front();
        q.pop();

        for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
        {
            const int &next = graph[current][i_edge];
            if (d[next] == -1)
            {
                q.push(next);
                d[next] = d[current] + 1;
                if (d[next] < mins[next])
                {
                    mins[next] = d[next];
                    argmins[next] = p;
                }
            }
        }
    }
}
vector<vector<int> > build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    vector<std::set<int> > undirected(n);
    vector<vector<int> > graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if edge not seen
        {
            undirected[i].insert(j);
            undirected[j].insert(i);
            graph[i].push_back(j);
            graph[j].push_back(i);
        }
    }
    return graph;
}
vector<int> maxmin_random_sp_unweighted(const vector<vector<int> >& graph, int n_pivots, int p0, int seed)
{
    int n = graph.size();

    vector<int> mins(n, std::numeric_limits<int>::max());
    vector<int> argmins(n, -1);

    // first pivot
    mins[p0] = 0;
    argmins[p0] = p0;
    maxmin_bfs_unweighted(graph, p0, mins, argmins);
    for (int i = 0; i < n; i++)
    {
        if (argmins[i] == -1)
            throw std::invalid_argument("graph has multiple connected components");
    }

    // remaining pivots
    // std::mt19937 rng(seed);
    // std::uniform_real_distribution<> uniform(0, 1);
    rk_state rstate;
    rk_seed(seed, &rstate);
    for (int ij = 1; ij < n_pivots; ij++)
    {
        // choose pivots with probability min
        int min_total = 0;
        for (int i = 0; i < n; i++)
        {
            min_total += mins[i];
        }
        int sample = rk_interval(min_total, &rstate);
        int cumul = 0;
        int argmax = -1;
        for (int i = 0; i < n; i++)
        {
            cumul += mins[i];
            if (cumul >= sample)
            {
                argmax = i;
                break;
            }
        }
        if (argmax == -1)
            throw std::invalid_argument("unweighted pivot sampling failed");

        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs_unweighted(graph, argmax, mins, argmins);
    }
    return argmins;
}
// is not actually a multi-source shortest path, because regions come for free with maxmin_random_sp
vector<constraint> MSSP_unweighted(const vector<vector<int> >& graph, const vector<int>& closest_pivots, int flag)
{
    int n = graph.size();

    // get pivots and their regions, but in sets
    std::map<int, std::set<int> > regions;
    std::map<int, std::map<int, constraint> > termsDict;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(closest_pivots[i]) == regions.end())
        {
            regions[closest_pivots[i]] = std::set<int>();
        }
        regions[closest_pivots[i]].insert(i);
    }

    std::map<int, std::set<int> >::iterator region;
    for (region=regions.begin(); region!=regions.end(); ++region)
    {
        // q contains next to visit
        std::queue<int> q;
        vector<int> d(n, -1);

        int p = region->first;
        q.push(p);
        d[p] = 0;

        // q2 contains visited vertices' distances for s calculation
        std::queue<int> q2;
        int s = 0;
        q2.push(0);

        while (!q.empty())
        {
            int current = q.front();
            q.pop();

            for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
            {
                const int &next = graph[current][i_edge];
                if (d[next] == -1)
                {
                    q.push(next);
                    d[next] = d[current] + 1;

                    // empty the second queue enough to calculate s
                    while (!q2.empty() && q2.front() <= d[next]/2)
                    {
                        q2.pop();
                        s += 1;
                    }
                    if (region->second.find(next) != region->second.end())
                    {
                        q2.push(d[next]);
                    }

                    int i = next;
                    if (i < p)
                    {
                        if (termsDict.find(i) == termsDict.end())
                            termsDict[i] = std::map<int, constraint>();
                        if (termsDict[i].find(p) == termsDict[i].end())
                            termsDict[i].insert(std::pair<int, constraint>(p, constraint(i, p, d[next])));

                        // termsDict[i].at(p).w_ij = s / ((double)d[next] * d[next]);
                        if(flag==2)termsDict[i].find(p)->second.wij = s / ((double)d[next] * d[next]);
                        else if(flag==3)termsDict[i].find(p)->second.wij = s / (double)d[next];
                    }
                    else
                    {
                        if (termsDict.find(p) == termsDict.end())
                            termsDict[p] = std::map<int, constraint>();
                        if (termsDict[p].find(i) == termsDict[p].end())
                            termsDict[p].insert(std::pair<int, constraint>(i, constraint(p, i, d[next])));

                        // termsDict[p].at(i).w_ji = s / ((double)d[next] * d[next]);
                        if(flag==2)termsDict[p].find(i)->second.wji = s / ((double)d[next] * d[next]);
                        else if(flag==3)termsDict[p].find(i)->second.wji = s / (double)d[next];
                    }
                }
            }
        }
    }
    // 1-stress
    if(flag==2){
        for (int i=0; i<n; i++)
        {
            for (unsigned i_edge=0; i_edge<graph[i].size(); i_edge++)
            {
                const int j = graph[i][i_edge];
                if (i < j)
                {
                    if (termsDict.find(i) == termsDict.end())
                        termsDict[i] = std::map<int, constraint>();
                    if (termsDict[i].find(j) == termsDict[i].end())
                        termsDict[i].insert(std::pair<int, constraint>(j, constraint(i, j, 1)));
                    else
                        // termsDict[i].at(j).d = 1;
                        termsDict[i].find(j)->second.d = 1;

                    // termsDict[i].at(j).w_ij = termsDict[i].at(j).w_ji = 1;
                    termsDict[i].find(j)->second.w = termsDict[i].find(j)->second.w = 1;
                }
            }
        }
    }
    vector<constraint> terms;
    std::map<int, std::map<int, constraint> >::iterator it;
    for (it=termsDict.begin(); it!=termsDict.end(); ++it)
    {
        std::map<int, constraint>::iterator jt;
        for (jt=it->second.begin(); jt!=it->second.end(); ++jt)
        {
            terms.push_back(jt->second);
        }
    }
    return terms;
}

void DFS1(int n, int m, const vector<vector<int>> &e, int st, vector<int> &kn, int neighbor) {
    bool visited[n];
    int layer[n];
    memset(layer,0,n);
    memset(visited,0,n);
    visited[st]=1;
    queue<int> q;
    q.push(st);
    layer[st]=0;
    kn.pop_back();
    while(!q.empty()){
        int u=q.front();
        q.pop();
        for(int i=0;i<e[u].size();i++){
            int v=e[u][i];
            if(!visited[v]) {
                layer[v]=layer[u]+1;
                visited[v]=1;
                if(layer[v]<=neighbor){
                    q.push(v);
                    kn.push_back(v);
                }
            }
        }
    }
}

vector<constraint> MSSP_unweighted_framework(const vector<vector<int> >& graph, const vector<int>& closest_pivots, int neighbor, vector<Link> links, vector<vector<double> > shortPat)
{
    int n = graph.size();
    int m=links.size();
    // get pivots and their regions, but in sets
    std::map<int, std::set<int> > regions;
    std::map<int, std::map<int, constraint> > termsDict;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(closest_pivots[i]) == regions.end())
        {
            regions[closest_pivots[i]] = std::set<int>();
        }
        regions[closest_pivots[i]].insert(i);
    }

    std::map<int, std::set<int> >::iterator region;
    for (region=regions.begin(); region!=regions.end(); ++region)
    {
        // q contains next to visit
        std::queue<int> q;
        vector<int> d(n, -1);

        int p = region->first;
        q.push(p);
        d[p] = 0;

        // q2 contains visited vertices' distances for s calculation
        std::queue<int> q2;
        int s = 0;
        q2.push(0);

        while (!q.empty())
        {
            int current = q.front();
            q.pop();

            for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
            {
                const int &next = graph[current][i_edge];
                if (d[next] == -1)
                {
                    q.push(next);
                    d[next] = d[current] + 1;

                    // empty the second queue enough to calculate s
                    while (!q2.empty() && q2.front() <= d[next]/2)
                    {
                        q2.pop();
                        s += 1;
                    }
                    if (region->second.find(next) != region->second.end())
                    {
                        q2.push(d[next]);
                    }

                    int i = next;
                    if (i < p)
                    {
                        if (termsDict.find(i) == termsDict.end())
                            termsDict[i] = std::map<int, constraint>();
                        if (termsDict[i].find(p) == termsDict[i].end())
                            termsDict[i].insert(std::pair<int, constraint>(p, constraint(i, p, d[next])));

                        // termsDict[i].at(p).w_ij = s / ((double)d[next] * d[next]);
                        termsDict[i].find(p)->second.wij = s;

                    }
                    else
                    {
                        if (termsDict.find(p) == termsDict.end())
                            termsDict[p] = std::map<int, constraint>();
                        if (termsDict[p].find(i) == termsDict[p].end())
                            termsDict[p].insert(std::pair<int, constraint>(i, constraint(p, i, d[next])));

                        // termsDict[p].at(i).w_ji = s / ((double)d[next] * d[next]);
                        termsDict[p].find(i)->second.wji = s ;

                    }
                }
            }
        }
    }
    vector<constraint> terms;

    // 1-stress
    if(neighbor==1){
        for (int i = 0; i < m; i++) {
            int src, tgt;
            src = links[i].source;
            tgt = links[i].target;
            if (src == tgt) continue;
            terms.push_back(constraint(src, tgt, 1, 1));
        }
    } else if(neighbor>=2){
        vector<Link> slink;
        int E1=0;

        for(int i=0;i<n;i++){
            vector<int> kn;
            kn.push_back(i);
            DFS1(n,m,graph,i,kn,neighbor);//�ĳ�n���ھ�
            //cout<<kn.size()<<endl;
            for(int j=0;j<kn.size();j++){
                //cout<<kn[j]<<" "<<i<<endl;
                if(kn[j]<i){
                    struct Link temp;
                    temp.source=i;
                    temp.target=kn[j];
                    slink.push_back(temp);
                }
            }
        }
        E1=slink.size();
        //cout<<E1<<" "<<m<<endl;
        for (int i = 0; i < slink.size(); i++) {
            int src, tgt;
            src = slink[i].source;
            tgt = slink[i].target;
            if (src == tgt) continue;

            terms.push_back(constraint(src, tgt, shortPat[src][tgt], 1));
        }
    }
    std::map<int, std::map<int, constraint> >::iterator it;
    for (it=termsDict.begin(); it!=termsDict.end(); ++it)
    {
        std::map<int, constraint>::iterator jt;
        for (jt=it->second.begin(); jt!=it->second.end(); ++jt)
        {
            terms.push_back(jt->second);
        }
    }

    //cout<<"neighbor"<<neighbor<<" "<<terms.size()<<endl;
    return terms;
}
void fisheryates_shuffle(vector<constraint> &terms, rk_state &rstate)
{
    int n = terms.size();
    //clock_t time_start = clock();
    for (unsigned i=n-1; i>=1; i--)
    {

        unsigned j = rk_interval(i, &rstate);
        constraint temp = terms[i];

        terms[i] = terms[j];
        terms[j] = temp;
        //std::swap(terms[i],terms[j]);
    }

    //clock_t time_end = clock();
    //cout << "rk time:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << endl;
}
void fisheryates_shuffle(vector<constraint_sgd> &terms, rk_state &rstate)
{
    int n = terms.size();
    //clock_t time_start = clock();
    for (unsigned i=n-1; i>=1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        constraint_sgd temp = terms[i];

        terms[i] = terms[j];
        terms[j] = temp;
        //std::swap(terms[i],terms[j]);
    }
}
void fisheryates_shuffle(vector<int> &terms, rk_state &rstate)
{
    int n = terms.size();
    // clock_t time_start = clock();
    for (unsigned i=n-1; i>=1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        int temp = terms[i];
        terms[i] = terms[j];
        terms[j] = temp;
        //std::swap(terms[i],terms[j]);
    }
//    clock_t time_end = clock();
//    cout << "rk time:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << endl;
}

void graph::solveDij() {

    vector<Edge> e[n];
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
        vector<float> dist(n);
        Dij(n,m,e,i,dist);
        vector<double> isp(n);
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
float Jaccard(set<int> set1,set<int> set2){
    set<int> inter_set;
    set<int> union_set;
    //sort(set1.begin(),set1.end());
    //sort(set2.begin(),set2.end());

    set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(), inserter(inter_set,inter_set.begin()));
    set_union(set1.begin(),set1.end(),set2.begin(),set2.end(), inserter(union_set,union_set.begin()));
    return (float)inter_set.size()/(float)union_set.size();
}

void DFS2(int n,int m,vector<Edge> e[MAX],int st,vector<int> &kn, int neighbor){
//void DFS_K(int n, int m, const vector<vector<int>> &e, int st, vector<int> &kn, int neighbor) {
    bool visited[n];
    int layer[n];
    memset(layer,0,n);
    memset(visited,0,n);
    visited[st]=1;
    queue<int> q;
    q.push(st);
    layer[st]=0;
    kn.pop_back();
    while(!q.empty()){
        int u=q.front();
        q.pop();
        for(int i=0;i<e[u].size();i++){
            int v=e[u][i].end;
            if(!visited[v]) {
                layer[v]=layer[u]+1;
                visited[v]=1;
                if(layer[v]<=neighbor){
                    q.push(v);
                    kn.push_back(v);
                }
            }
        }

    }
}

void mult(vector<float>& b_k, vector<vector<float>>& B, int k) {
    vector<float> b_k1(k);
    for (int i = 0; i < k; i++) {
        b_k1[i] = b_k[i];
        b_k[i] = 0;
    }
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            b_k[i] += B[i][j] * b_k1[j];
        }
    }
}
float norm(vector<float>& b_k, int k) {
    float a = 0;
    for (int i = 0; i < k; i++) {
        a += (b_k[i] * b_k[i]);
    }
    return a;
}
void PivotMDS(vector<Node> &nodes,int  k,int N,vector<vector<double>> &shortPat,vector<int> &pivot){
    int iteration = 100;
    vector<vector<float>> C(N, vector<float>(k));
    vector<vector<int>> dkj(N, vector<int>(k));
    //select k pivots
    //rand
    //set ShortestPath
    //calculate B
    //calculate CC^T
    //calculate first second

    //Init dik
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            dkj[i][j] = shortPat[i][pivot[j]];

        }
    }
    float sum = 0;
    vector<float> sumn(N);
    vector<float> sumk(k);
    for (int i = 0; i < N; i++) {
        sumn[i] = 0;
        // sumk[i] = 0;
    }
    for (int i = 0; i < k; i++) {
        // sumn[i] = 0;
        sumk[i] = 0;
    }
    //calculate Cij
    for (int i = 0; i < N; i++) {
        //sumk[i] = 0;
        for (int j = 0; j < k; j++) {
            //dkj[i][j]=1;//shortPat[pivot[i]][j]
            sum += (dkj[i][j] * dkj[i][j]);
            sumn[i] += (dkj[i][j] * dkj[i][j]);
        }
        sumn[i] /= k;
    }
    sum /= (N * k);

    for (int j = 0; j < k; j++) {
        //sumn[j] = 0;
        for (int i = 0; i < N; i++) {
            //dkj[i][j]=1;//shortPat[pivot[i]][j]
            //sum+=(dkj[i][j]*dkj[i][j]);
            sumk[j] += (dkj[i][j] * dkj[i][j]);
        }
        sumk[j] /= N;
    }
    //calculate C
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j] = -0.5f * ((dkj[i][j] * dkj[i][j]) - sumk[j] - sumn[i] + sum);
        }
    }
    vector<vector<float>> B(k, vector<float>(k));
    // B.resize(N,vector<float>(N));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            for (int _ = 0; _ < N; _++) {
                B[i][j] += C[_][i] * C[_][j];
            }
        }
    }
    vector<float> b_k(k);
    for (int i = 0; i < k; i++) {
        if (i % 2 == 0) {
            b_k[i] = 1;
        }
        else {
            b_k[i] = -1;
        }
    }
    //?????????????
    for (int i = 0; i < iteration; i++) {
        mult(b_k, B, k);
        //float norm1=norm(b_k,k);
        float a = 0;
        for (int i = 0; i < k; i++) {
            a += (b_k[i] * b_k[i]);
        }
        for (int i = 0; i < k; i++) {
            b_k[i] /=sqrtf(a);
        }
    }
    //?????
    float normV = norm(b_k, k);
    float a = 0;
    for (int i = 0; i < k; i++) {
        a += (b_k[i] * b_k[i]);
    }
    for (int i = 0; i < k; i++) {
        b_k[i] /= sqrtf(a);
    }
    //???????
    float lamdab = 0;
    vector<float> temp(k);
    for (int i = 0; i < k; i++) {
        temp[i] = b_k[i];
    }
    mult(temp, B, k);

    for (int i = 0; i < k; i++) {
        lamdab += temp[i] * b_k[i];
        // cout << "temp" << i << "=   b_k="<<b_k[i] << " ";
    }
    vector<vector<float>> A(k, vector<float>(k));
    vector<float> a_k(k);
    float lamdaa = 0;
    // float normV = norm(b_k, k);

    //??A
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A[i][j] = B[i][j] - lamdab * b_k[i] * b_k[j] / normV;
        }
    }
    //?????A_k
    for (int i = 0; i < k; i++) {
        if (i % 2 == 0) {
            a_k[i] = 1;
        }
        else {
            a_k[i] = -1;
        }
    }

    //??��?????????
    for (int i = 0; i < iteration; i++) {
        mult(a_k, A, k);
        //float norm1=norm(b_k,k);
        float a2 = 0;
        for (int i = 0; i < k; i++) {
            a2 += (a_k[i] * a_k[i]);
        }
        for (int i = 0; i < k; i++) {
            a_k[i] /= sqrtf(a2);
        }
    }
    //?????
    float a2 = 0;
    for (int i = 0; i < k; i++) {
        a2 += (a_k[i] * a_k[i]);
    }
    //norm
    //ak
    for (int i = 0; i < k; i++) {
        a_k[i] /= sqrtf(a2);
    }

    // vector<float> x(N);
    // vector<float> y(N);

    for (int i = 0; i < N; i++) {
        nodes[i].x = 0;
        nodes[i].y = 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            nodes[i].x  += C[i][j] * b_k[j];
            nodes[i].y += C[i][j] * a_k[j];
        }
    }
    /*for (int i = 0; i < N; i++) {
        cout<<nodes[i].x<<" "<<nodes[i].y<<" ";
    }*/


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
                cout<<"pivot"<<endl;
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
            if(src<tgt){
                //constraints.push_back(constraint(src,tgt,sp[src][tgt],w,pa));
                int temp=src;
                src=tgt;
                tgt=temp;
            }
            auto it = set_sgd.find(constraint_sgd(src,tgt));
            if(it!=set_sgd.end()){
                it->force_sgd.push_back(pa);
            }else{
                vector<force> f;
                f.push_back(pa);
                set_sgd.insert(constraint_sgd(src,tgt,sp[src][tgt],w,f));
            }

            /*if(pivots.size()){
                cout<<"pivot"<<endl;
                constraints.push_back(constraint(tgt,src,sp[src][tgt],w,pa));//np
            }*/

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
                //constraints.push_back(constraint(i,j,sp[i][j],w,pa));
                auto it = set_sgd.find(constraint_sgd(i,j));
                if(it!=set_sgd.end()){
                    it->force_sgd.push_back(pa);
                }else{
                    vector<force> f;
                    f.push_back(pa);
                    set_sgd.insert(constraint_sgd(i,j,sp[i][j],w,f));
                }
            }

        }
    }

}
void graph::append_range(vector<Link> range,para pa){
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
void graph::append_range(vector<Link> range,force pa){
    for(int i=0;i<range.size();i++){
        int src=range[i].source,tgt=range[i].target;
        if(src!=tgt){
            float w;
            if(sp[src][tgt]!=0){
                w= 1.0f / (sp[src][tgt] * sp[src][tgt]);
            }else{
                w = 0.0;
            }
            if(src<tgt){
                int temp=src;
                src=tgt;
                tgt=temp;
            }
            auto it = set_sgd.find(constraint_sgd(src,tgt));
            if(it!=set_sgd.end()){
                it->force_sgd.push_back(pa);
            }else{
                vector<force> f;
                f.push_back(pa);
                set_sgd.insert(constraint_sgd(src,tgt,sp[src][tgt],w,f));
            }
        }
    }
}
void graph::append_Np_range(para pa) {

    //consSGD.insert()
    vector<constraint> constraintnp;
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
    vector<vector<int> > g1 = build_graph_unweighted(n, m, I, J);
    vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 42);
    constraintnp= MSSP_unweighted_framework(g1, closest_pivots,2,links,sp);
    for(int i=0;i<closest_pivots.size();i++){
        pivots[closest_pivots[i]]=1;
    }

    set<int> pset;
    for(int i=0;i<closest_pivots.size();i++){
        pset.insert(closest_pivots[i]);
    }
    vector<int> pivot;
    /*randperm(p,n,pivot);
    sort(pivot.begin(),pivot.end());*/
    for(set<int>::iterator it=pset.begin();it!=pset.end();it++){
        pivot.push_back(*it);
        // cout<<*it<<" ";
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
    vector<constraint> constraintnp;
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
    vector<vector<int> > g1 = build_graph_unweighted(n, m, I, J);
    vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 42);
    constraintnp= MSSP_unweighted_framework(g1, closest_pivots,2,links,sp);
    for(int i=0;i<closest_pivots.size();i++){
        pivots[closest_pivots[i]]=1;
    }

    set<int> pset;
    for(int i=0;i<closest_pivots.size();i++){
        pset.insert(closest_pivots[i]);
    }
    vector<int> pivot;
    /*randperm(p,n,pivot);
    sort(pivot.begin(),pivot.end());*/
    for(set<int>::iterator it=pset.begin();it!=pset.end();it++){
        pivot.push_back(*it);
        // cout<<*it<<" ";
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
                auto it = set_sgd.find(constraint_sgd(src,tgt));
                if(it!=set_sgd.end()){
                    it->force_sgd.push_back(pa);
                }else{
                    vector<force> f;
                    f.push_back(pa);
                    set_sgd.insert(constraint_sgd(src,tgt,sp[src][tgt],w,f));
                }
            }
        }
    }
    /*for(int i=0;i<constraintnp.size();i++){
        consSGD.insert(constraintnp[i]);
    }*/

}
void graph::append_S_range(int neighbor,para pa) {
    //consSGD.insert()
    vector<Edge> e[n];
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
        vector<int> kn;
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
    vector<Edge> e[n];
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
        vector<int> kn;
        kn.push_back(i);
        DFS2(n,m,e,i,kn,neighbor);
        for(int j=0;j<kn.size();j++){
            if(kn[j]<i){
                float w = 1.0f / (sp[i][kn[j]] * sp[i][kn[j]]);
                auto it = set_sgd.find(constraint_sgd(i,kn[j]));
                if(it!=set_sgd.end()){
                    it->force_sgd.push_back(pa);
                }else{
                    vector<force> f;
                    f.push_back(pa);
                    set_sgd.insert(constraint_sgd(i,kn[j],sp[i][kn[j]],w,f));
                }
            }
        }
    }
}
bool graph::power_law_graph() {
    vector<int> mask(m+1,0);
    vector<int> count(n,0);
    int maxD=0;
    bool isPow=false;
    for(int i=0;i<m;i++){
        count[links[i].source]++;
        count[links[i].target]++;
    }
    for(int i=0;i<this->n;i++){
        mask[count[i]]++;
        maxD=max(maxD,mask[count[i]]);
    }
    if(mask[1]>0.8*maxD&&mask[1]>0.3*m) {isPow= true;cout<<"graph is low power graph,set p=1.8"<<endl;}
    return isPow;
}

void graph::preSolveSGD(double eps,int t_max,int seed,float eta_max){
    float w_min = 100000000;
    float w_max = 0;
    for(int k=0;k<constraints.size();k++){
        float tempw=constraints[k].w;
        w_min = min(tempw, w_min);
        w_max = max(tempw, w_max);
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
        vector<vector<int> > g1 = build_graph_unweighted(n, m, I, J);
        vector<int> closest_pivots = maxmin_random_sp_unweighted(g1, p, 0, 48);
        set<int> pset;
        for(int i=0;i<closest_pivots.size();i++){
            pset.insert(closest_pivots[i]);
        }
        vector<int> pivot;
        //randperm(p,n,pivot);
        //sort(pivot.begin(),pivot.end());
        for(set<int>::iterator it=pset.begin();it!=pset.end();it++){
            pivot.push_back(*it);
            // cout<<*it<<" ";
        }
        /*for(int i=0;i<pivot.size();i++){
            cout<<pivot[i]<<",";
        }*/
        cout<<"compute pivotMDS"<<endl;
        PivotMDS(nodes,pivot.size(),n,sp,pivot);
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
        //cout<<"schedule"<<schedule[i]<<endl;
    }
}
void graph::preSolveSGD2(double eps,int t_max,int seed,float eta_max){

    constraints_sgd.assign(set_sgd.begin(),set_sgd.end());
    float w_min = 100000000;
    float w_max = 0;
    for(int k=0;k<constraints_sgd.size();k++){
        float tempw=constraints_sgd[k].w;
        w_min = min(tempw, w_min);
        w_max = max(tempw, w_max);
    }
    rk_seed(seed, &rstate);
    fisheryates_shuffle(constraints_sgd, rstate);
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
    //cout << "time ran:" << 1000 * (time_end - time_start) / (double)CLOCKS_PER_SEC << " " <<iter<< endl;
    double Delta_max=0;
    // cout<<"11111"<<endl;
//#pragma omp parallel for schedule(dynamic)
//cout<<nodes[0].x<<","<<nodes[0].y<<endl;
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
            cout<<"q="<<q<<"p="<<p<<"a="<<a<<"r="<<r<<"ka="<<ka<<"kr="<<kr<<endl;
        }*/
        const float &pm=t.pa.pmds;
        const int &i = t.i, &j = t.j;
        const double &w_ij = t.w;
        const double &d_ij = t.d;
        if(d_ij==0) continue;
        alpha = schedule[iter];
        //cout<<con<<" get schedule"<<endl;
        //if(alpha>1.0f)alpha=1.0f;
        float mvx = (nodes[i].x - nodes[j].x);
        float mvy = (nodes[i].y - nodes[j].y);
        float dist = sqrtf(mvx * mvx + mvy * mvy);
        //cout<<con<<" get dist"<<endl;
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
                //cout<<"q="<<q<<endl;
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
                //cout<<"a="<<a<<endl;
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
                //cout<<"p="<<p<<endl;
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

            //if(iter==0) cout<<"no pivot"<<endl;
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
    // cout<<"11111"<<endl;
//#pragma omp parallel for schedule(dynamic)
//cout<<nodes[0].x<<","<<nodes[0].y<<endl;
    for (int con = 0; con < constraints_sgd.size(); con++) {
        float rx = 0;
        float ry = 0;

        const constraint_sgd &t = constraints_sgd[con];
        const int &i = t.i, &j = t.j;
        const double &w_ij = t.w;
        const double &d_ij = t.d;
        alpha = schedule[iter];
        /*if(iter==0){
            cout<<"i="<<i<<"j="<<j<<"w="<<w_ij<<"d="<<d_ij<<endl;
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
             cout<<"k="<<k<<"a="<<a<<"b="<<b<<endl;*/
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
                //cout<<"a="<<a<<endl;
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
    vector<Edge> e[n];
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
        vector<int> kn;
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
    //cout<<force[0]<<" "<<force[1]<<" ";
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
        //cout<<force[0]<<" "<<force[1]<<" ";
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
void Quad::add(float x, float y, int idx) {
    size++;
    float tempx = mass[0], tempy = mass[1];
    bool sameflag = false;
    float mult1 = (float)(size - 1.0f) / size;
    float mult2 = (float)1.0f / size;
    mass[0] = mass[0] * mult1 + mult2 * x;
    mass[1] = mass[1] * mult1 + mult2 * y;
    if (is_leaf & !have_node) {
        have_node = true;
        nodeindex = idx;
        return;
    }

    if (is_leaf) {
        if ((x - mass[0]) * (x - mass[0]) + (y - mass[1]) * (y - mass[1]) > 1e-12) {
            children[0] = new Quad(center[0] - 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[1] = new Quad(center[0] + 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[2] = new Quad(center[0] - 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            children[3] = new Quad(center[0] + 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            is_leaf = false;
            if (mult1 <= 0.0f) {
                printf("In here must more than one node and mult1>0.0f");
                exit(1);
            }
            int index = 0;
            if (tempx > center[0])index += 1;
            if (tempy > center[1])index += 2;
            children[index]->add(tempx, tempy, nodeindex);

        }
        else {
            children[0] = new Quad(center[0] - 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[1] = new Quad(center[0] + 0.25f * w, center[1] - 0.25f * h, 0.5f * w, 0.5f * h);
            children[2] = new Quad(center[0] - 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            children[3] = new Quad(center[0] + 0.25f * w, center[1] + 0.25f * h, 0.5f * w, 0.5f * h);
            is_leaf = false;
            if (mult1 <= 0.0f) {
                printf("In here must more than one node and mult1>0.0f");
                exit(1);
            }
            int index = 0;
            if (tempx > center[0])index += 1;
            if (tempy > center[1])index += 2;
            children[index]->add(tempx, tempy, nodeindex);
            sameflag = true;
        }
    }

    if (!is_leaf) {
        int index = 0;
        if (x > center[0])index += 1;
        if (y > center[1])index += 2;
        index = (index + sameflag) % 4;
        children[index]->add(x, y, idx);
    }
}
Quad::Quad(vector<Node>& nodes, int N) {
    float maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9;
    for (int i = 0; i < N; i++) {
        maxx = nodes[i].x > maxx ? nodes[i].x : maxx;
        minx = nodes[i].x < minx ? nodes[i].x : minx;
        maxy = nodes[i].y > maxy ? nodes[i].y : maxy;
        miny = nodes[i].y < miny ? nodes[i].y : miny;
    }
    is_leaf = true;
    have_node = false;
    nodeindex = -1;
    size = 0;
    mass[0] = mass[1] = .0;
    center[0] = (maxx + minx) / 2.0f;
    center[1] = (maxy + miny) / 2.0f;

    w = (maxx - minx) > (maxy - miny) ? (maxx - minx) : (maxy - miny);
    w = 1.01f * w;
    h = w;
    children = (Quad**)malloc(4 * sizeof(Quad*));
    children[0] = children[1] = children[2] = children[3] = NULL;
    for (int i = 0; i < N; i++) {
        add(nodes[i].x, nodes[i].y, i);
    }

}

Quad::Quad(float cx, float cy, float width, float height) {
    is_leaf = true;
    have_node = false;
    size = 0;
    mass[0] = mass[1] = .0;
    center[0] = cx, center[1] = cy;
    w = width;
    h = height;
    nodeindex = -1;
    children = (Quad**)malloc(4 * sizeof(Quad*));
    children[0] = children[1] = children[2] = children[3] = NULL;
}

Quad::~Quad() {
    for (int i = 0; i < 4; i++) {
        if (children[i] != NULL)
            delete children[i];
    }
    free(children);
}

void QuadForce(Quad* qt, float& fx, float& fy, float posx, float posy, float theta, float alpha, int nodeindex,float dp,int layer) {
    //fprintf(stderr,"Hello QuadForce\n");
    //bool debugFlag = false;
    float eps=0;
    if(dp>1.0){eps=0.1;}else{eps=0.01;}
    if (qt->size <= 0.0)return;
    float mvx = (posx - qt->mass[0]);
    float mvy = (posy - qt->mass[1]);
    float dist = mvx * mvx + mvy * mvy;
    // if(debugFlag)fprintf(stderr,"%.16f %f %f\n",dist,mvx,mvy);
    if (qt->is_leaf && dist <= 1e-4) {
        if (qt->size > 1) {
            printf("No qtnode have more than one node");
            exit(1);
        }
        //if (debugFlag)fprintf(stderr, "%d %d\n", qt->nodeindex, nodeindex);
        if (qt->nodeindex != nodeindex) {
            if (mvx == 0.0) { mvx = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvx * mvx; }
            if (mvy == 0.0) { mvy = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvy * mvy; }
            if (dist < 1e-12) {
                dist = sqrtf(0.01f * dist);
            }
            // if(debugFlag)fprintf(stderr,"dist,%f %f %f %f %f %f %f %f %f\n",qt->mass[0],qt->mass[1],qt->w,qt->size,dist,mvx,mvy,fx,fy);
            //********* force define  *********//
            dist=sqrtf(dist) + eps;
            fx -= alpha * mvx * (qt->size - 1) / pow(dist,dp+1);
            fy -= alpha * mvy * (qt->size - 1) / pow(dist,dp+1);
        }
    }
    else if (qt->is_leaf || dist / (qt->w * qt->w) > theta * theta) {

        if (mvx == 0.0) { mvx = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvx * mvx; }
        if (mvy == 0.0) { mvy = 0.001f * (rand() / RAND_MAX - 0.5f), dist += mvy * mvy; }
        if (dist <= 1e-5) {
            dist = sqrtf(0.01f * dist) + eps;
        }
        // if(debugFlag)fprintf(stderr,"theta,%f %f %f %f %f %f %f %f %f\n",qt->mass[0],qt->mass[1],qt->w,qt->size,dist,mvx,mvy,fx,fy);
        //********* force define  *********//
        dist=sqrtf(dist) + eps;
        fx -= alpha * mvx * qt->size /pow(dist,dp+1);
        fy -= alpha * mvy * qt->size /pow(dist,dp+1);
    }
    else {
        layer++;
       if(layer>12){return;}
        QuadForce(qt->children[0], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
        QuadForce(qt->children[1], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
        QuadForce(qt->children[2], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
        QuadForce(qt->children[3], fx, fy, posx, posy, theta, alpha, nodeindex,dp,layer);
    }
}
float testNP3(vector<Node> &nodes,vector<Link> &links,vector<vector<double>> &sp,int neigh){
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
    //cout << "NP" <<neigh<<"="<<NP / number << endl;
    return NP / number;
}

/*
void graph::outToFile(string filename,string methodname) {

    string str1= ".\\" + methodname;
    //system(str1.c_str());
    if (access(str1.c_str(),6)==-1)
    {
        mkdir(str1.c_str());
    }
    string str2= ".\\" + methodname+"\\"+ filename.substr(0, filename.length() - 4);
    if (access(str2.c_str(),6)==-1)
    {
        mkdir(str2.c_str());
    }
    string filename1 = methodname+"\\"+ filename.substr(0, filename.length() - 4)+"\\"+"nodes.csv";
    freopen((char*) filename1.data(),"w",stdout);
    //freopen(str,"w",stdout);
    cout<<"positionx"<<","<<"positiony"<<endl;

    for(int i=0;i<this->n;i++){
        cout<<nodes[i].x<<","<<nodes[i].y<<endl;

    }

    //freopen("/dev/tty", "w", stdout);
    freopen("CON", "r", stdin);
    freopen("CON", "w", stdout);
}

void graph::RelationOutToFile(string filename,string methodname) {

    string str1= ".\\" + methodname;
    //system(str1.c_str());
    if (access(str1.c_str(),6)==-1)
    {
        mkdir(str1.c_str());
    }
    string str2= ".\\" + methodname+"\\"+ filename.substr(0, filename.length() - 4);
    if (access(str2.c_str(),6)==-1)
    {
        mkdir(str2.c_str());
    }
    string filename1 = methodname+"\\"+ filename.substr(0, filename.length() - 4)+"\\"+"links.csv";
    freopen((char*) filename1.data(),"w",stdout);
    cout<<"source"<<","<<"target"<<endl;
    for(int i=0;i<this->m;i++){
        cout<<links[i].source<<","<<links[i].target<<endl;
    }
    //freopen("/dev/tty", "w", stdout);
    freopen("CON", "r", stdin);
    freopen("CON", "w", stdout);
}
*/
void graph::initPosition(string filename) {
    ifstream input(filename);
    string s;
    //cout<<"read data"<<endl;
    if(nodes.size()!=n){
        nodes.resize(n);
    }
    int index=0;

    while(getline(input, s)){
        string s1,s2;
        istringstream is(s);
        is>>s1>>s2;
        nodes[index].x= str2float(s1);
        nodes[index].y= str2float(s2);
        index++;
    }
}
void graph::drawSVG(string filename, string method) {
    string filename1 = method+"_" +filename.substr(0, filename.length() - 4)+"_layout.svg";
    float maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9;
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
        nodes[i].x*=1000;
        nodes[i].x+=100;
        nodes[i].y/=lenth;
        nodes[i].y*=1000+100;
        nodes[i].y+=100;
    }
    fprintf(stderr,"new drawSVG\n");
    freopen((char*) filename1.data(),"w",stdout);
    cout<<"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n"
          "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
          "\n"
          "<svg width=\"1200\" height=\"1200\" version=\"1.1\"\n"
          "xmlns=\"http://www.w3.org/2000/svg\">"<<endl;
    
    for(int i=0;i<this->m;i++){
        //cout<<links[i].source<<","<<links[i].target<<endl;
       cout<<"<line x1=\""<<nodes[links[i].source].x<<"\" y1=\""<<nodes[links[i].source].y<<"\" x2=\""<<nodes[links[i].target].x<<"\" y2=\""<<nodes[links[i].target].y<<"\"\n"
             "style=\"stroke:rgb(99,99,99);stroke-width:2\"/>"<<endl;
    }
    for(int i=0;i<this->n;i++){
        //cout<<links[i].source<<","<<links[i].target<<endl;
        cout<<"<circle cx=\""<<nodes[i].x<<"\" cy=\""<<nodes[i].y<<"\" r=\"3\" stroke=\"green\"\n"
              "stroke-width=\"1.5\" fill=\"blue\"/>"<<endl;
    }
    cout<<"</svg>"<<endl;
    
    //freopen("CON", "r", stdin);
    //freopen("CON", "w", stdout);
    freopen("/dev/tty", "w", stdout);
}