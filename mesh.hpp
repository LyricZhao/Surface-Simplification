# ifndef __MESH_H__
# define __MESH_H__

# define MAX_LENGTH 65536ul

# include <cmath>
# include <fstream>
# include <iostream>
# include <set>
# include <string>
# include <vector>

char buffer[MAX_LENGTH];

inline double sqr(double x) { return x * x; }

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
        inline void add_face(const Face &face) {
            faces.insert(face);
            return;
        }
    };

    int tot;
    std:: vector<Vertex> vertexes; /* All indexes start from zero. */

    void get_p(const Face &face, double *p);
    void add_kp(const Face &face, double *kp);
    int face_count();
    int get_father(int v);
    std:: string read_number(std:: string line, int &pos);
    void add_vertex(std:: string line);
    void add_face(std:: string line);
    void calculate_Q();
    void select_pairs();
    void aggregation();
    void refresh_memory();

public:
    Mesh();
    ~Mesh();

    void simplify(double ratio);
    void load_from_file(std:: string path);
    void write_into_file(std:: string path);
};

Mesh:: Mesh() { }

Mesh:: ~Mesh() { }

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

void Mesh:: calculate_Q() {
    for(auto &v: vertexes) {
        memset(v.q, 0, sizeof(int) * 16);
        for(auto &face: v.faces)
            add_kp(face, v.q);
    }
    return;
}

void Mesh:: simplify(double ratio) {
    calculate_Q();
    select_pairs();
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