# ifndef __MESH_H__
# define __MESH_H__

# define MAX_LENGTH 65536ul

# include <cmath>
# include <fstream>
# include <iostream>
# include <map>
# include <queue>
# include <set>
# include <string>
# include <vector>

char buffer[MAX_LENGTH];
const double eps = 1e-8;

inline double sqr(double x) {
    return x * x;
}

class Mesh {
private:
    struct Face {
        int dim[3];

        Face(int *v) {
            dim[0] = v[0], dim[1] = v[1], dim[2] = v[2];
            return;
        }

        inline void sort() {
            int min_d = std:: min(std:: min(dim[0], dim[1]), dim[2]);
            if(min_d == dim[0])
                return;
            else if(min_d == dim[1]) /* 012 -> 120 */
                std:: swap(dim[0], dim[1]), std:: swap(dim[1], dim[2]);
            else if(min_d == dim[2]) /* 012 -> 201 */
                std:: swap(dim[0], dim[1]), std:: swap(dim[0], dim[2]);
            return;
        }
        inline friend std:: ostream& operator << (std:: ostream &os, const Face &face) {
            os << "f " << face.dim[0] << " " << face.dim[1] << " " << face.dim[2];
            return os;
        }
        inline bool operator < (const Face &b) const {
            return dim[0] != b.dim[0] ? (dim[0] < b.dim[0]) :
                (dim[1] != b.dim[1] ? (dim[1] < b.dim[1]) :
                (dim[2] < b.dim[2]));
        }
    };

    struct Vertex {
        int father;
        double dim[3], q[16];
        std:: set<Face> faces;

        inline friend std:: ostream& operator << (std:: ostream &os, const Vertex &v) {
            os << "v " << v.dim[0] << " " << v.dim[1] << " " << v.dim[2];
            return os;
        }
        inline double distance(const Vertex &v) {
            return sqrt(sqr(dim[0] - v.dim[0]) + sqr(dim[1] - v.dim[1]) + sqr(dim[2] - v.dim[2]));
        }
        inline void add_face(const Face &face) {
            faces.insert(face);
            return;
        }
    };

    struct Pair {
        int v0, v1, time;
        double cost, dim[3];

        inline void sort() { if(v0 > v1) std:: swap(v0, v1); return; }
        inline bool operator < (const Pair &b) const {
            return cost < b.cost;
        }
        inline unsigned long long index_hash() const {
            return (((unsigned long long)v0) << 32) | ((unsigned long long)v1);
        }
    };

    struct Heap {
        /* The heap will only maintain the newest cost. */
        int stamp;
        std:: map<unsigned long long, int> times;
        std:: priority_queue<Pair> queue;
        
        Heap() { stamp = 0; return; }
        inline void clear() {
            stamp = 0, times.clear();
            while(queue.size())
                queue.pop();
            return;
        }
        inline bool count(unsigned long long hash) {
            return times.count(hash) ? (times[hash] > 0) : false;
        }
        inline void refresh() {
            while(queue.size() > 0) {
                const Pair &pair = queue.top();
                int time = times[pair.index_hash()];
                if(time == pair.time) break;
                queue.pop();
            }
            return;
        }
        inline Pair top() {
            refresh();
            return queue.top();
        }
        inline void pop() {
            refresh();
            queue.pop();
            return;
        }
        inline void insert(Pair &pair) {
            pair.time = ++ stamp;
            times[pair.index_hash()] = stamp;
            return;
        }
        inline void del(const Pair &pair) {
            times[pair.index_hash()] = -1;
            return;
        }
    };

    int tot, face_tot;
    std:: vector<Vertex> vertexes; /* All indexes start from zero. */
    std:: vector<std:: set<int> > edges;
    Heap heap;

    void add_pair(int v0, int v1);
    void select_dfs(int v, int p, int depth, double t, std:: set<int> &searched);
    double get_cost(double *dim, double *q);
    void get_new_v(double *q, double *dim, double *dim_a, double *dim_b);
    void get_p(const Face &face, double *p);
    void add_kp(const Face &face, double *kp);
    int face_count();
    int get_father(int v);
    std:: string read_number(std:: string line, int &pos);
    void add_vertex(std:: string line);
    void add_face(std:: string line);
    void calculate_Q();
    void select_pairs(double t);
    void aggregation();
    void refresh_memory();

public:
    Mesh();
    ~Mesh();

    void simplify(double ratio, double t);
    void load_from_file(std:: string path);
    void write_into_file(std:: string path);
};

Mesh:: Mesh() { tot = face_tot = 0; return; }

Mesh:: ~Mesh() { }

int Mesh:: face_count() {
    return face_tot;
}

int Mesh:: get_father(int v) {
    int &father = vertexes[v].father;
    return father == v ? v : father = get_father(father);
}

std:: string Mesh:: read_number(std:: string line, int &pos) {
    auto is_number = [] (char c) -> bool {
        return ('0' <= c && c <= '9') || c == '.';
    };
    int start = pos;
    while(start < line.length() && !is_number(line[start])) ++ start;
    pos = start + 1;
    while(pos < line.length() && is_number(line[pos])) ++ pos;
    return line.substr(start, pos - start);
}

void Mesh:: add_vertex(std:: string line) {
    Vertex vertex;
    vertex.father = tot ++;
    for(int i = 0, pos = 0; i < 3; ++ i)
        vertex.dim[i] = atof(read_number(line, pos).c_str());
    return;
}

void Mesh:: add_face(std:: string line) {
    ++ face_tot;
    int v[3];
    bool read_twice = (line.find('/') != line.npos);
    for(int i = 0, pos = 0; i < 3; ++ i) {
        v[i] = atoi(read_number(line, pos).c_str()) - 1;
        if(read_twice) read_number(line, pos);
    }
    Face face(v); face.sort();
    for(int i = 0; i < 3; ++ i)
        vertexes[v[i]].add_face(face);
    return;
}

double Mesh:: get_cost(double *dim, double *q) {
    double q[16], t[4], nd[4];
    for(int i = 0; i < 4; ++ i) t[i] = 0;
    nd[0] = dim[0], nd[1] = dim[1], nd[2] = dim[2], nd[3] = 1.;
    for(int i = 0; i < 4; ++ i) for(int j = 0; j < 4; ++ j)
        t[i] += nd[j] * q[j * 4 + i];
    double cost = 0;
    for(int i = 0; i < 4; ++ i)
        cost += t[i] * nd[i];
    return cost;
}

void Mesh:: get_p(const Face &face, double *p) {
    Vertex &a = vertexes[face.dim[0]], &b = vertexes[face.dim[1]], &c = vertexes[face.dim[2]];
    p[0] = (b.dim[1] - a.dim[1]) * (c.dim[2] - a.dim[2]) - (c.dim[1] - a.dim[1]) * (b.dim[2] - a.dim[2]);
    p[1] = (b.dim[2] - a.dim[2]) * (c.dim[0] - a.dim[0]) - (c.dim[2] - a.dim[2]) * (b.dim[0] - a.dim[0]);
    p[2] = (b.dim[0] - a.dim[0]) * (c.dim[1] - a.dim[1]) - (c.dim[0] - a.dim[0]) * (b.dim[1] - a.dim[1]);
    double div = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
    p[0] /= div, p[1] /= div, p[2] /= div;
    p[3] = -p[0] * a.dim[0] - p[1] * a.dim[1] - p[2] * a.dim[2];
    return;
}

void Mesh:: add_kp(const Face &face, double *kp) {
    double p[4];
    get_p(face, p);
    for(int i = 0; i < 4; ++ i) for(int j = 0; j < 4; ++ j)
        kp[i * 4 + j] += p[i] * p[j];
    return;
}

void Mesh:: get_new_v(double *q, double *dim, double *dim_a, double *dim_b) {
    double d = q[0] * q[5] * q[10] - q[0] * q[6] * q[6] - q[1] * q[1] * q[10] + q[1] * q[6] * q[2] + q[2] * q[1] * q[6] - q[2] * q[5] * q[2];
    if(d > eps) {
        double x = q[3] * q[5] * q[10] - q[3] * q[6] * q[6] - q[7] * q[1] * q[10] + q[7] * q[6] * q[2] + q[11] * q[1] * q[6] - q[11] * q[5] * q[2];
        double y = q[0] * q[7] * q[10] - q[0] * q[11] * q[6] - q[1] * q[3] * q[10] + q[1] * q[11] * q[2] + q[2] * q[3] * q[6] - q[2] * q[7] * q[2];
        double z = q[0] * q[5] * q[11] - q[0] * q[6] * q[7] - q[1] * q[1] * q[11] + q[1] * q[6] * q[3] + q[2] * q[1] * q[7] - q[2] * q[5] * q[3];
        dim[0] = -x / d, dim[1] = -y / d, dim[2] = -z / d;
    } else {
        dim[0] = (dim_a[0] + dim_b[0]) / 2., dim[1] = (dim_a[1] + dim_b[1]) / 2., dim[2] = (dim_a[2] + dim_b[2]) / 2.;
    }
    return;
}

void Mesh:: add_pair(int v0, int v1) {
    if(v0 > v1) std:: swap(v0, v1);
    Pair pair; pair.v0 = v0, pair.v1 = v1;
    Vertex &a = vertexes[v0], &b = vertexes[v1];
    double q[16];
    for(int i = 0; i < 16; ++ i) q[i] = a.q[i] + b.q[i];
    get_new_v(q, pair.dim, a.dim, b.dim);
    pair.cost = get_cost(pair.dim, q);
    heap.insert(pair);
    return;
}

void Mesh:: select_dfs(int v, int p, int depth, double t, std:: set<int> &searched) {
    if((depth > 1 && vertexes[v].distance(vertexes[p].distance) > t) || searched.count(p) > 0) return;
    searched.insert(p);
    if(depth >= 0 && v != p)
        add_pair(v, p);
    for(auto next: edges[v])
        select_dfs(v, next, depth + 1, t, searched);
    return;
}

void Mesh:: select_pairs(double t) {
    for(int i = 0; i < vertexes.size(); ++ i) {
        for(auto face: vertexes[i].faces) {
            for(int j = 0; j < 3; ++ j) if(face.dim[j] != i)
                edges[i].insert(face.dim[j]);
        }
    }
    for(int i = 0; i < vertexes.size(); ++ i) {
        std:: set<int> searched;
        select_dfs(i, i, 0, t, searched);
    }
    return;
}

void Mesh:: calculate_Q() {
    for(auto &v: vertexes) {
        memset(v.q, 0, sizeof(int) * 16);
        for(auto &face: v.faces)
            add_kp(face, v.q);
    }
    return;
}

void Mesh:: aggregation() {
    Pair pair = heap.top(); heap.pop();

    return;
}

void Mesh:: refresh_memory() {
    heap.clear();
    edges.clear();
    return;
}

void Mesh:: simplify(double ratio, double t) {
    calculate_Q();
    select_pairs(t);
    int target = face_count() * ratio;
    while(face_count() > target)
        aggregation();
    refresh_memory();
    return;
}

void Mesh:: load_from_file(std:: string path) {
    std:: ifstream input(path);
    while(input.getline(buffer, MAX_LENGTH)) {
        std:: string line = buffer;
        if(!line.length()) continue;
        if(line[0] == 'v') add_vertex(line);
        if(line[0] == 'f') add_face(line);
    }
    return;
}

void Mesh:: write_into_file(std:: string path) {
    int *index = (int*) std:: malloc(sizeof(int) * tot), cnt = 0;
    std:: ofstream output(path);
    for(int i = 0; i < vertexes.size(); ++ i) {
        if(get_father(i) != i) continue;
        index[i] = ++ cnt;
        output << vertexes[i] << std:: endl;
    }
    for(int i = 0; i < vertexes.size(); ++ i) {
        if(get_father(i) != i) continue;
        for(auto &face: vertexes[i].faces) {
            if(face.dim[0] != i) continue;
            output << face << std:: endl;
        }
    }
    std:: free(index);
    return;
}

# endif