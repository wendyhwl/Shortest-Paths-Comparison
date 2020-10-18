#include <ctime>
#include <cmath>
#include <queue>
#include <vector>
#include <cctype>
#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <algorithm>
using namespace std;

const double eps = 1e-7;
const double INF = DBL_MAX / 3;
using VVI = vector<vector<int>>;

// define node
struct Node {
    int id;
    double x;
    double y;
};

// graph class
class Graph {
  public:

    // Graph(input file)
    Graph(const char *file) : _v(0) {
        char buf[128];
        int id;
        double x, y;

        // open input file & read each line
        FILE *fp = fopen(file, "r");
        while (fgets(buf, 128, fp)) {
            if (strlen(buf) == 1) break;
            sscanf(buf, "%d : %lf, %lf", &id, &x, &y);
            // assume node id is from 1 to n
            // thus _nodes[i].id == i + 1
            _nodes.push_back({id, x, y});
            ++_v;
        }
        
        // resize container so it contains _v element
        _g.resize(_v);
        while (fgets(buf, 128, fp)) {
            int i = 0;
            int u = 0;
            for (; buf[i]!=':'; ++i) {
                if (isdigit(buf[i])) u = u * 10 + buf[i] - '0';
            }
            ++i;
            char *cp = strtok(buf+i, ", ");
            while (cp) {
                if (cp && cp[0]) {
                    int v = atoi(cp);
                    if (v > 0) _g[u-1].push_back(v-1);
                }
                cp = strtok(NULL, ",");
            }
        }
        fclose(fp);
    }
    
    // getV()
    int getV() const { return _v; }

    // &getVs(int u)
    const vector<int> &getVs(int u) const { return _g[u]; }

    // distance(int u, int v)
    // caculate distance between vertices
    double distance(int u, int v) const {
        const Node &q1 = _nodes[u];
        const Node &q2 = _nodes[v];
        double dlat = 2 * M_PI * (q2.x - q1.x) / 360;
        double mlat = 2 * M_PI * (q1.x + q2.x) / 2 / 360;
        double dlon = 2 * M_PI * (q2.y - q1.y) / 360;
        double x = cos(mlat) * dlon;
        double y = dlat;
        return sqrt(x*x+y*y) * 6371009;
    }

  private:
    int _v;
    VVI _g;
    vector<Node> _nodes;
};

// equal
bool equal(double x, double y) { return fabs(x-y) < eps; }

// dijkstra
pair<int, double> dijkstra(const Graph &g, int src, int tgt) {
    int n = g.getV();

    // initiate dis storing total distance
    vector<double> dis(n, INF+1);

    // create priority queue q
    priority_queue<pair<double, int>> q;

    // distance to source set as 0, push to p
    dis[src] = 0;
    q.push({0, src});
    int vn = 0;

    // while q is not empty, extract min from q
    while (!q.empty()) {
        auto &pr = q.top();
        double w = -pr.first;
        int u = pr.second;
        q.pop();

        // initiate visited node vn
        ++vn;

        if (!equal(dis[u], w)) continue;
        if (u == tgt) break;

        // for each unvisited neighbour, update its distance
        for (auto v: g.getVs(u)) {
            double d = g.distance(u, v);

            // update dis iff its original dis is larger
            if (dis[v] > w + d) {
                dis[v] = w + d;
                q.push({-dis[v], v});
            }
        }
    }

    // return number of visited vertices and distance between src and tgt
    return {vn, dis[tgt]};
}

// astar
pair<int, double> astar(const Graph &g, int src, int tgt) {
    int n = g.getV();

    // initiate dis storing total distance
    vector<double> dis(n, INF+1);

    // { {dis[u]+to_t[u], dis[u]}, u}
    // creating priority queue q
    priority_queue<pair<double, pair<double, int>>> q;

    // distance to source set as 0, push to q
    dis[src] = 0;
    q.push({-g.distance(src, tgt), {0, src}});

    // initiate visited node vn
    int vn = 0;

    // while q is not empty
    while (!q.empty()) {

        // current node pr = node with least dis value
        auto &pr = q.top();
        double w = -pr.second.first;
        int u = pr.second.second;

        // extract current node pr
        q.pop();

        // initiate visited node vn
        ++vn;

        if (!equal(dis[u], w)) continue;

        // done if the current node pr is the goal
        if (u == tgt) break;

        for (auto v: g.getVs(u)) {
            double d = g.distance(u, v);

            // for testing
            // printf("distance: %lf\n", d);

            // update dis iff its original dis is larger
            if (dis[v] > w + d) {
                dis[v] = w + d;
                q.push({-dis[v]-g.distance(v, tgt), {-dis[v], v}});
            }
        }
    }

    // return number of visited vertices and distance between src and tgt
    return {vn, dis[tgt]};
}

// generate_marks
// given Graph g, Source src
// for each vertex v in g, calcuate the shortest distance from src to v
// stored in dis, that is dis[v] = shortest-path-from-src-to-v
void generate_marks(const Graph &g, int src, vector<double> &dis) {
    // using dijkstra algorithm
    int n = g.getV();
    dis.resize(n, INF);
    priority_queue<pair<double, int>> q;

    dis[src] = 0;
    q.push({0, src});
    while (!q.empty()) {
        auto &pr = q.top();
        double w = -pr.first;
        int u = pr.second;
        q.pop();
        if (!equal(dis[u], w)) continue;
        for (auto v: g.getVs(u)) {
            double d = g.distance(u, v);
            if (dis[v] > w + d) {
                dis[v] = w + d;
                q.push({-dis[v], v});
            }
        }
    }
}

// landmark
pair<int, double> landmark(const Graph &g, int src, int tgt) {
    int n = g.getV();
    int step = n / 4;
    // create 3 landmarks
    vector<int> marks = {step, 2*step, 3*step};
    vector<vector<double>> mark_distances(marks.size());
    // generate distance for each landmark
    for (int i=0,len=(int)marks.size(); i<len; ++i) {
        generate_marks(g, marks[i], mark_distances[i]);
    }

    // helper function to calcuate estimate distance
    // = max(|dist[u]-dist[t]|)
    auto to_t = [&mark_distances, tgt](int u) -> double {
        double mx = -INF;
        for (auto &dist: mark_distances) {
            double d = fabs(dist[u] - dist[tgt]);
            mx = max(mx, d);
        }
        return mx;
    };

    vector<double> dis(n, INF+1);
    // { {dis[u]+max(|dist[u]-dist[t]|), dis[u]}, u}
    priority_queue<pair<double, pair<double, int>>> q;

    // same procedure like A*
    // the only difference is to_t
    dis[src] = 0;
    q.push({-to_t(src), {0, src}});
    int vn = 0;
    while (!q.empty()) {
        auto &pr = q.top();
        double w = -pr.second.first;
        int u = pr.second.second;
        q.pop();
        ++vn;
        if (!equal(dis[u], w)) continue;
        if (u == tgt) break;
        for (auto v: g.getVs(u)) {
            double d = g.distance(u, v);

            // for testing
            // printf("distance: %lf\n", d);
            if (dis[v] > w + d) {
                dis[v] = w + d;
                q.push({-dis[v]-to_t(v), {-dis[v], v}});
            }
        }
    }

    // return number of visited vertices and distance between src and tgt
    return {vn, dis[tgt]};
}

// main
int main(int argc, char **argv) {
    Graph g("graph1000.txt");
    int n = g.getV();
    srand(time(0));
    int num = 20;
    int nd = 0, na = 0, nl = 0;

    // generate 20 random sets and printing their results
    for (int i=0; i<num; ++i) {
        // get two random number from 1-1000
        int q1 = rand() % n;
        int q2 = rand() % n;

        printf("QUERY %d: %d -> %d\n", i+1, q1+1, q2+1);
        // calculating distance & number of nodes visited
        const auto &_d = dijkstra(g, q1, q2);
        const auto &_a = astar(g, q1, q2);
        const auto &_l = landmark(g, q1, q2);

        // store number of nodes visited in variables _nd, _na, _nl
        int _nd = _d.first;
        int _na = _a.first;
        int _nl = _l.first;

        // printing results for all three methods
        if (_d.second >= INF) printf("1. Dijkstra: %-5d (INF)\n", _nd);
        else printf("1. Dijkstra: %-5d (%lf)\n", _nd, _d.second);
        if (_a.second >= INF) printf("2. Astar:    %-5d (INF)\n", _na);
        else printf("2. Astar:    %-5d (%lf)\n", _na, _a.second);
        if (_l.second >= INF) printf("3. Landmark: %-5d (INF)\n", _nl);
        else printf("3. Landmark: %-5d (%lf)\n", _nl, _l.second);
        puts("");

        // store result in nd, na, nl for calculating average
        nd += _nd;
        na += _na;
        nl += _nl;
    }

    // printing average number of nodes visited
    printf("Average\n");
    printf("1. Dijkstra %lf\n", double(nd)/num);
    printf("2. Astar    %lf\n", double(na)/num);
    printf("3. Landmark %lf\n", double(nl)/num);
    return 0;

} // end of code