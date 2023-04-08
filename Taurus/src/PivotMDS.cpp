#include "PivotMDS.h"
void Dij(int V,int E,std::vector<Edge> e[MAX],int st,std::vector<float> &dist){
    std::priority_queue<std::Pair,std::vector<std::Pair>,std::greater<std::Pair>> q;
    int visited[V];
    memset(visited,0,sizeof(visited));
    for(int i=0;i<V;i++) dist[i]=MAX;
    dist[st]=0;
    q.push(std::make_pair(0,st));
    while(!q.empty()){
        std::Pair t=q.top();
        q.pop();
        if(visited[t.second]) continue;
        visited[t.second]=1;

        for(int i=0;i<e[t.second].size();i++){
            int son=e[t.second][i].end;
            if(dist[son]>dist[t.second]+e[t.second][i].w){
                dist[son]=dist[t.second]+e[t.second][i].w;
                if(!visited[son])
                    q.push(std::make_pair(dist[son],son));
            }
        }
    }
}
void maxmin_bfs_unweighted(const std::vector<std::vector<int> >& graph, const int p, std::vector<int>& mins, std::vector<int>& argmins)
{
    int n = graph.size();
    std::queue<int> q;
    std::vector<int> d(n, -1);

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
std::vector<std::vector<int> > build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    std::vector<std::set<int>> undirected(n);
    std::vector<std::vector<int> > graph(n);

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
std::vector<int> maxmin_random_sp_unweighted(const std::vector<std::vector<int> >& graph, int n_pivots, int p0, int seed)
{
    int n = graph.size();

    std::vector<int> mins(n, std::numeric_limits<int>::max());
    std::vector<int> argmins(n, -1);

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
std::vector<constraint> MSSP_unweighted(const std::vector<std::vector<int> >& graph, const std::vector<int>& closest_pivots, int flag)
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
        std::vector<int> d(n, -1);

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
    std::vector<constraint> terms;
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

void DFS1(int n, int m, const std::vector<std::vector<int>> &e, int st, std::vector<int> &kn, int neighbor) {
    bool visited[n];
    int layer[n];
    memset(layer,0,n);
    memset(visited,0,n);
    visited[st]=1;
    std::queue<int> q;
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

std::vector<constraint> MSSP_unweighted_framework(const std::vector<std::vector<int> >& graph, const std::vector<int>& closest_pivots, int neighbor, std::vector<Link> links, std::vector<std::vector<double> > shortPat)
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
        std::vector<int> d(n, -1);

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
    std::vector<constraint> terms;

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
        std::vector<Link> slink;
        int E1=0;

        for(int i=0;i<n;i++){
            std::vector<int> kn;
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

void BFS(int V,int E,std::vector<Edge> e[MAX],int st,std::vector<int> &ind,std::vector<int>& count_com){
    //int visited[V];
    int *visited=new int[V];
    memset(visited,0,sizeof(visited));
    std::queue<int> q1;
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
float testNP3(std::vector<Node> &nodes,std::vector<Link> &links,std::vector<std::vector<double>> &sp,int neigh){
    double NP = 0;
    int leafnum = 0;
    int number=nodes.size();
    int countt=0;
    for (int i = 0; i < number; i++) {
        std::vector<int> dijknn;
        std::vector<double> dknn;
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
float Jaccard(std::set<int> set1,std::set<int> set2){
    std::set<int> inter_set;
    std::set<int> union_set;
    //sort(set1.begin(),set1.end());
    //sort(set2.begin(),set2.end());

    std::set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(), inserter(inter_set,inter_set.begin()));
    std::set_union(set1.begin(),set1.end(),set2.begin(),set2.end(), inserter(union_set,union_set.begin()));
    return (float)inter_set.size()/(float)union_set.size();
}

void DFS2(int n,int m,std::vector<Edge> e[MAX],int st,std::vector<int> &kn, int neighbor){
//void DFS_K(int n, int m, const vector<vector<int>> &e, int st, vector<int> &kn, int neighbor) {
    bool visited[n];
    int layer[n];
    memset(layer,0,n);
    memset(visited,0,n);
    visited[st]=1;
    std::queue<int> q;
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

void mult(std::vector<float>& b_k, std::vector<std::vector<float>>& B, int k) {
    std::vector<float> b_k1(k);
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
float norm(std::vector<float>& b_k, int k) {
    float a = 0;
    for (int i = 0; i < k; i++) {
        a += (b_k[i] * b_k[i]);
    }
    return a;
}
void PivotMDS(std::vector<Node> &nodes,int  k,int N,std::vector<std::vector<double>> &shortPat,std::vector<int> &pivot){
    int iteration = 100;
    std::vector<std::vector<float>> C(N, std::vector<float>(k));
    std::vector<std::vector<int>> dkj(N, std::vector<int>(k));
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
    std::vector<float> sumn(N);
    std::vector<float> sumk(k);
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
    std::vector<std::vector<float>> B(k, std::vector<float>(k));
    // B.resize(N,vector<float>(N));

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            for (int _ = 0; _ < N; _++) {
                B[i][j] += C[_][i] * C[_][j];
            }
        }
    }
    std::vector<float> b_k(k);
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
    std::vector<float> temp(k);
    for (int i = 0; i < k; i++) {
        temp[i] = b_k[i];
    }
    mult(temp, B, k);

    for (int i = 0; i < k; i++) {
        lamdab += temp[i] * b_k[i];
        // cout << "temp" << i << "=   b_k="<<b_k[i] << " ";
    }
    std::vector<std::vector<float>> A(k, std::vector<float>(k));
    std::vector<float> a_k(k);
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