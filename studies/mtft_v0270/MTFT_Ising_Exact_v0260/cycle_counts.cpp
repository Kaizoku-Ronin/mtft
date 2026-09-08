// Independent exact route: enumerate every binary cycle-space element.
// Input: vertex_count edge_count, followed by edge endpoint pairs.
// Output: JSON array of even-subgraph counts by number of occupied edges.
// Uses neither the Ising Pfaffian code nor polynomial variable elimination.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

struct Bits { uint64_t lo=0, hi=0; };
Bits operator^(Bits a, Bits b) { return {a.lo^b.lo,a.hi^b.hi}; }
Bits bit(int e) { return e<64 ? Bits{uint64_t(1)<<e,0} : Bits{0,uint64_t(1)<<(e-64)}; }
struct DSU {
    std::vector<int> p;
    explicit DSU(int n):p(n) { std::iota(p.begin(),p.end(),0); }
    int root(int v) { return p[v]==v ? v : p[v]=root(p[v]); }
    bool join(int u,int v) { u=root(u);v=root(v);if(u==v)return false;p[u]=v;return true; }
};

int main() {
    int n,E;
    if(!(std::cin>>n>>E) || n<=0 || E<0 || E>128) return 2;
    auto started=std::chrono::steady_clock::now();
    std::vector<std::pair<int,int>> edges(E);
    std::vector<std::vector<std::pair<int,int>>> tree(n);
    std::vector<int> chords;
    DSU dsu(n);
    for(int e=0;e<E;++e) {
        auto &[u,v]=edges[e];std::cin>>u>>v;
        if(u<0 || v<0 || u>=n || v>=n) return 2;
        if(dsu.join(u,v)) { tree[u].push_back({v,e});tree[v].push_back({u,e}); }
        else chords.push_back(e);
    }
    if(chords.size()>29) { std::cerr<<"Cycle rank exceeds verification budget 29\n";return 3; }
    std::vector<Bits> cycles;
    for(int e:chords) {
        auto [u,v]=edges[e];
        std::vector<int> parent(n,-1), pe(n,-1);
        std::queue<int> q;q.push(u);parent[u]=u;
        while(!q.empty() && parent[v]<0) {
            int a=q.front();q.pop();
            for(auto [b,f]:tree[a]) if(parent[b]<0) { parent[b]=a;pe[b]=f;q.push(b); }
        }
        if(parent[v]<0) throw std::runtime_error("missing tree path");
        Bits C=bit(e);
        for(int a=v;a!=u;a=parent[a]) C=C^bit(pe[a]);
        cycles.push_back(C);
    }
    uint64_t count=uint64_t(1)<<cycles.size();
    std::vector<uint64_t> histogram(E+1,0);histogram[0]=1;
    Bits state;
    for(uint64_t i=1;i<count;++i) {
        state=state^cycles[__builtin_ctzll(i)];
        ++histogram[__builtin_popcountll(state.lo)+__builtin_popcountll(state.hi)];
    }
    std::cout<<"[";
    for(int e=0;e<=E;++e) std::cout<<(e ? "," : "")<<histogram[e];
    std::cout<<"]\n";
    std::chrono::duration<double> elapsed=std::chrono::steady_clock::now()-started;
    std::cerr<<"Enumerated "<<count<<" even subgraphs (rank "<<cycles.size()<<") in "<<elapsed.count()<<" seconds\n";
}
