# ifndef __MESH_H__
# define __MESH_H__

# define MAX_LENGTH 65536ul

# define DFS_T_SEARCH
# define FLIP_COST
# define SHARP_COST

# include <algorithm>
# include <cassert>
# include <cstring>
# include <cmath>
# include <fstream>
# include <iostream>
# include <map>
# include <queue>
# include <set>
# include <string>
# include <vector>

char buffer[MAX_LENGTH];
const double eps = 1e-12;

inline double sqr(double x) {
    return x * x;
}

class Mesh {
private:
    struct Face {
        /* Note: the order can not be ignored */
        int dim[3];

        Face(int *v) {
            dim[0] = v[0], dim[1] = v[1], dim[2] = v[2];
            sort();
            return;
        }

        Face(int v0, int v1, int v2) {
            dim[0] = v0, dim[1] = v1, dim[2] = v2;
            sort();
            return;
        }

        inline void replace(int v0, int v1) {
            if(dim[0] == v0) dim[0] = v1;
            else if(dim[1] == v0) dim[1] = v1;
            else if(dim[2] == v0) dim[2] = v1;
            else assert(0);
            return;
        }
        inline int find_diff(int v0) {
            if(dim[0] != v0) return dim[0];
            if(dim[1] != v0) return dim[1];
            if(dim[2] != v0) return dim[2];
            assert(0);
        }
        inline int find_diff(int v0, int v1) {
            if(dim[0] != v0 && dim[0] != v1) return dim[0];
            if(dim[1] != v0 && dim[1] != v1) return dim[1];
            if(dim[2] != v0 && dim[2] != v1) return dim[2];
            assert(0);
        }
        inline bool contain(int v) const {
            return dim[0] == v || dim[1] == v || dim[2] == v;
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
        friend bool operator < (const Face &a, const Face &b) {
            return a.dim[0] != b.dim[0] ? (a.dim[0] < b.dim[0]) :
                (a.dim[1] != b.dim[1] ? (a.dim[1] < b.dim[1]) :
                (a.dim[2] < b.dim[2]));
        }
    };

    struct Vertex {
        int father;
        double dim[3], q[16];
        std:: set<Face> faces;

        inline void replace(double *nd, double *qb) {
            dim[0] = nd[0], dim[1] = nd[1], dim[2] = nd[2];
            for(int i = 0; i < 16; ++ i) q[i] += qb[i];
            return;
        }
        inline std:: vector<Face> find_another(int id, int v, bool option=false) {
            std:: vector<Face> another;
            for(auto &face: faces) {
                if(face.contain(v) == option) continue;
                another.push_back(face);
            }
            return another;
        }
        inline bool count_face(const Face &face) {
            return faces.count(face) > 0;
        }
        inline void remove_face(const Face &face) {
            faces.erase(face);
            return;
        }
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

        Pair(int _v0, int _v1) {
            v0 = _v0, v1 = _v1;
            sort();
            return;
        }

        inline void sort() { if(v0 > v1) std:: swap(v0, v1); return; }
        inline bool operator < (const Pair &b) const {
            return cost > b.cost;
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
        std:: vector<std:: set<int> > in_queue_pairs;
        
        Heap() { stamp = 0; return; }
        inline void resize_queue_pairs(int n) {
            in_queue_pairs.resize(n);
            return;
        }
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
            queue.push(pair);
            in_queue_pairs[pair.v0].insert(pair.v1);
            in_queue_pairs[pair.v1].insert(pair.v0);
            return;
        }
        inline void del(const Pair &pair) {
            times[pair.index_hash()] = -1;
            in_queue_pairs[pair.v0].erase(pair.v1);
            in_queue_pairs[pair.v1].erase(pair.v0);
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
    double cos_distance(double *d1, double *d2);
    void get_norm(double *dim, double *a_dim, double *b_dim, double *c_dim);
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
        return ('0' <= c && c <= '9') || c == '.' || c == '-' || c == 'e';
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
    vertexes.push_back(vertex);
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
    Face face(v);
    for(int i = 0; i < 3; ++ i)
        vertexes[v[i]].add_face(face);
    return;
}

double Mesh:: get_cost(double *dim, double *q) {
    double t[4], nd[4];
    for(int i = 0; i < 4; ++ i) t[i] = 0;
    nd[0] = dim[0], nd[1] = dim[1], nd[2] = dim[2], nd[3] = 1.;
    for(int i = 0; i < 4; ++ i) for(int j = 0; j < 4; ++ j)
        t[i] += nd[j] * q[j * 4 + i];
    double cost = 0;
    for(int i = 0; i < 4; ++ i)
        cost += t[i] * nd[i];
    return cost;
}

void Mesh:: get_norm(double *dim, double *a_dim, double *b_dim, double *c_dim) {
    dim[0] = (b_dim[1] - a_dim[1]) * (c_dim[2] - a_dim[2]) - (c_dim[1] - a_dim[1]) * (b_dim[2] - a_dim[2]);
    dim[1] = (b_dim[2] - a_dim[2]) * (c_dim[0] - a_dim[0]) - (c_dim[2] - a_dim[2]) * (b_dim[0] - a_dim[0]);
    dim[2] = (b_dim[0] - a_dim[0]) * (c_dim[1] - a_dim[1]) - (c_dim[0] - a_dim[0]) * (b_dim[1] - a_dim[1]);
    return;
}

void Mesh:: get_p(const Face &face, double *p) {
    double *a_dim = vertexes[face.dim[0]].dim, *b_dim = vertexes[face.dim[1]].dim, *c_dim = vertexes[face.dim[2]].dim;
    p[0] = (b_dim[1] - a_dim[1]) * (c_dim[2] - a_dim[2]) - (c_dim[1] - a_dim[1]) * (b_dim[2] - a_dim[2]);
    p[1] = (b_dim[2] - a_dim[2]) * (c_dim[0] - a_dim[0]) - (c_dim[2] - a_dim[2]) * (b_dim[0] - a_dim[0]);
    p[2] = (b_dim[0] - a_dim[0]) * (c_dim[1] - a_dim[1]) - (c_dim[0] - a_dim[0]) * (b_dim[1] - a_dim[1]);
    double div = sqrt(sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
    p[0] /= div, p[1] /= div, p[2] /= div;
    p[3] = -p[0] * a_dim[0] - p[1] * a_dim[1] - p[2] * a_dim[2];
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
    if(fabs(d) > eps) {
        double x = q[3] * q[5] * q[10] - q[3] * q[6] * q[6] - q[7] * q[1] * q[10] + q[7] * q[6] * q[2] + q[11] * q[1] * q[6] - q[11] * q[5] * q[2];
        double y = q[0] * q[7] * q[10] - q[0] * q[11] * q[6] - q[1] * q[3] * q[10] + q[1] * q[11] * q[2] + q[2] * q[3] * q[6] - q[2] * q[7] * q[2];
        double z = q[0] * q[5] * q[11] - q[0] * q[6] * q[7] - q[1] * q[1] * q[11] + q[1] * q[6] * q[3] + q[2] * q[1] * q[7] - q[2] * q[5] * q[3];
        dim[0] = -x / d, dim[1] = -y / d, dim[2] = -z / d;
    } else {
        dim[0] = (dim_a[0] + dim_b[0]) / 2., dim[1] = (dim_a[1] + dim_b[1]) / 2., dim[2] = (dim_a[2] + dim_b[2]) / 2.;
    }
    return;
}

double Mesh:: cos_distance(double *d1, double *d2) {
    return (d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]) / (sqrt(sqr(d1[0]) + sqr(d1[1]) + sqr(d1[2])) * sqrt(sqr(d2[0]) + sqr(d2[1]) + sqr(d2[2])));
}

void Mesh:: add_pair(int v0, int v1) {
    if(v0 > v1) std:: swap(v0, v1);
    Pair pair(v0, v1);
    double *a_dim = vertexes[v0].dim, *b_dim = vertexes[v1].dim;
    double *a_q = vertexes[v0].q, *b_q = vertexes[v1].q;
    double q[16];
    for(int i = 0; i < 16; ++ i) q[i] = a_q[i] + b_q[i];
    get_new_v(q, pair.dim, a_dim, b_dim);
    pair.cost = get_cost(pair.dim, q);

    double na_dim[3], nb_dim[3], nc_dim[3], norm_0[3], norm_1[3];
    /* Flip case */
    # ifdef FLIP_COST
    std:: vector<Face> changed_faces[2];
    changed_faces[0] = vertexes[v0].find_another(v0, v1, true);
    changed_faces[1] = vertexes[v1].find_another(v1, v0, true);
    double avg_cos_d = 0; int num = 0;
    for(int c = 0; c < 2; ++ c) for(auto &face: changed_faces[c]) {
        memcpy(na_dim, vertexes[face.dim[0]].dim, sizeof(double) * 3);
        memcpy(nb_dim, vertexes[face.dim[1]].dim, sizeof(double) * 3);
        memcpy(nc_dim, vertexes[face.dim[2]].dim, sizeof(double) * 3);
        get_norm(norm_0, na_dim, nb_dim, nc_dim);
        if(face.dim[0] == v0 || face.dim[0] == v1) memcpy(na_dim, pair.dim, sizeof(double) * 3);
        if(face.dim[1] == v0 || face.dim[1] == v1) memcpy(nb_dim, pair.dim, sizeof(double) * 3);
        if(face.dim[2] == v0 || face.dim[2] == v1) memcpy(nc_dim, pair.dim, sizeof(double) * 3);
        get_norm(norm_1, na_dim, nb_dim, nc_dim);
        avg_cos_d += cos_distance(norm_0, norm_1);
        ++ num;
    }
    if(num) {
        avg_cos_d /= num; avg_cos_d += 1.1;
        pair.cost /= std:: min(1., avg_cos_d / 1.3);
    }
    # endif

    /* Sharp feature cost */
    # ifdef SHARP_COST
    std:: vector<Face> common_faces = vertexes[v0].find_another(v0, v1);
    if(common_faces.size() == 2) {
        int a, b, c;
        a = common_faces[0].dim[0], b = common_faces[0].dim[1], c = common_faces[0].dim[2];
        get_norm(norm_0, vertexes[a].dim, vertexes[b].dim, vertexes[c].dim);
        a = common_faces[1].dim[0], b = common_faces[1].dim[1], c = common_faces[1].dim[2];
        get_norm(norm_1, vertexes[a].dim, vertexes[b].dim, vertexes[c].dim);
        pair.cost /= std:: min(1., cos_distance(norm_0, norm_1) + 1.1);
    }
    # endif

    heap.insert(pair);
    return;
}

/* Potential bug: dfs is invalid when v0 and v1 are unconnected, and the distance is not accurate. */
void Mesh:: select_dfs(int v, int p, int depth, double t, std:: set<int> &searched) {
    double dist = vertexes[v].distance(vertexes[p]);
    if((depth > 1 && dist > t) || searched.count(p) > 0) return;
    searched.insert(p);
    if(depth >= 0 && v < p)
        add_pair(v, p);
    for(auto next: edges[p])
        select_dfs(v, next, depth + 1, t, searched);
    return;
}

void Mesh:: select_pairs(double t) {
    std:: cout << "Selecting pairs with t = " << t << " ... " << std:: flush;
    edges.resize(vertexes.size());
    heap.resize_queue_pairs(vertexes.size());
    for(int i = 0; i < vertexes.size(); ++ i) {
        for(auto face: vertexes[i].faces) {
            for(int j = 0; j < 3; ++ j) if(face.dim[j] != i)
                edges[i].insert(face.dim[j]);
        }
    }
    if(t < eps) {
        for(int i = 0; i < vertexes.size(); ++ i) for(int j: edges[i]) if(i < j)
            add_pair(i, j);
    } else {
        # ifdef DFS_T_SEARCH
            for(int i = 0; i < vertexes.size(); ++ i) {
                std:: set<int> searched;
                select_dfs(i, i, 0, t, searched);
            }
        # else
            std:: vector<int> index;
            index.resize(vertexes.size());
            for(int i = 0; i < index.size(); ++ i) index[i] = i;
            auto cmp = [this] (int i, int j) -> bool { return vertexes[i].dim[0] < vertexes[j].dim[0]; };
            std:: sort(index.begin(), index.end(), cmp);
            for(int i = 0; i < vertexes.size(); ++ i) for(int j: edges[i]) if(i < j)
                add_pair(i, j);
            for(int i = 0; i < vertexes.size(); ++ i) {
                int a = index[i], b;
                for(int j = i + 1; j < vertexes.size(); ++ j) {
                    b = index[j];
                    if(vertexes[a].distance(vertexes[b]) < t || !edges[a].count(b))
                        add_pair(a, b);
                    if(vertexes[b].dim[0] - vertexes[a].dim[0] > t) break;
                }
            }
        # endif
    }
    std:: cout << "ok !" << std:: endl;
    return;
}

void Mesh:: calculate_Q() {
    std:: cout << "Calculating Q matrix ... " << std:: flush;
    for(auto &v: vertexes) {
        memset(v.q, 0, sizeof(int) * 16);
        for(auto &face: v.faces)
            add_kp(face, v.q);
    }
    std:: cout << "ok !" << std:: endl;
    return;
}

void Mesh:: aggregation() {
    Pair pair = heap.top(); heap.pop();

    /* Merge set */
    int v0 = pair.v0, v1 = pair.v1;
    assert((get_father(v0) == v0) && (get_father(v1) == v1));
    vertexes[v1].father = v0;

    /* Remove the faces containing v0 and v1, no need to remove v1 */
    std:: vector<Face> another = vertexes[v0].find_another(v0, v1);
    for(auto face: another) {
        int v = face.find_diff(v0, v1);
        vertexes[v].remove_face(face);
        vertexes[v0].remove_face(face);
        -- face_tot;
    }

    /* Connect the connected points of v1 to v0, no need to remove v1 */
    another = vertexes[v1].find_another(v1, v0, true);
    for(auto original_face: another) {
        int a = original_face.find_diff(v1), b = original_face.find_diff(a, v1);
        Face new_face(original_face); new_face.replace(v1, v0);
        vertexes[a].remove_face(original_face), vertexes[b].remove_face(original_face);
        if(!vertexes[v0].count_face(new_face))
            vertexes[a].add_face(new_face), vertexes[b].add_face(new_face), vertexes[v0].add_face(new_face);
        else
            -- face_tot;
    }

    /* Replace v0 with the new node */
    vertexes[v0].replace(pair.dim, vertexes[v1].q);

    /* Remove all the pairs containing v1 in the queue, and connect them to v0 */
    std:: set<int> in_queue_pairs_0 = heap.in_queue_pairs[v0]; /* Note: here must be copying */
    std:: set<int> in_queue_pairs_1 = heap.in_queue_pairs[v1]; /* Note: here must be copying */
    for(auto v: in_queue_pairs_1) {
        heap.del(Pair(v, v1));
        if(!in_queue_pairs_0.count(v) && v != v0)
            add_pair(v, v0);
    }

    /* Update the value containing v0 in the queue. */
    for(auto v: in_queue_pairs_0) if(v != v1)
        add_pair(v, v0);
    
    /* TODO: considering new pair in the condition of t */
    return;
}

void Mesh:: simplify(double ratio, double t) {
    std:: cout << "Starting to simlify ... ok !" << std:: endl;
    calculate_Q();
    select_pairs(t);

    std:: cout << "Doing aggregation ... " << std:: flush;
    int target = face_count() * ratio;
    while(face_count() > target)
        aggregation();
    std:: cout << "ok!" << std:: endl;
    return;
}

void Mesh:: load_from_file(std:: string path) {
    std:: cout << "Loading from obj file " + path + " ... " << std:: flush;
    std:: ifstream input(path);
    while(input.getline(buffer, MAX_LENGTH)) {
        std:: string line = buffer;
        if(!line.length()) continue;
        if(line[0] == 'v') add_vertex(line);
        if(line[0] == 'f') add_face(line);
    }
    std:: cout << "ok !" << std:: endl;
    return;
}

void Mesh:: write_into_file(std:: string path) {
    std:: cout << "Writing into obj file " + path + " ... " << std:: flush;
    int *index = (int*) std:: malloc(sizeof(int) * tot), cnt = 0;
    std:: ofstream output(path);
    for(int i = 0; i < vertexes.size(); ++ i) {
        if(get_father(i) != i) continue;
        index[i] = ++ cnt;
        output << vertexes[i] << std:: endl;
    }
    double error_sum = 0; int vertex_cnt = 0, face_cnt = 0;
    for(int i = 0; i < vertexes.size(); ++ i) {
        if(get_father(i) != i) continue;
        ++ vertex_cnt;
        error_sum += get_cost(vertexes[i].dim, vertexes[i].q);
        for(auto &face: vertexes[i].faces) {
            if(face.dim[0] != i) continue;
            output << "f " << index[face.dim[0]] << " " << index[face.dim[1]] << " " << index[face.dim[2]] << std:: endl;
            ++ face_cnt;
        }
    }
    std:: free(index);
    std:: cout << "ok !" << std:: endl;
    std:: cout << "Vertex: " << vertex_cnt << std:: endl;
    std:: cout << "Face: " << face_cnt << std:: endl;
    std:: cout << "Avg error: " << error_sum / vertex_cnt << std:: endl;
    return;
}

# endif