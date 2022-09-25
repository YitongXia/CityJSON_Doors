#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <queue>
#include <cmath>

#include "DCEL.hpp"
#include "include/json.hpp"

using namespace std;
using json = nlohmann::json;

//class Point{
//public:
//    double x;
//    double y;
//    double z;
//
//    Point(double x1,double y1,double z1): x(x1),y(y1),z(z1)
//    {}
//};

void select_single_building(string &filename, string &building_num)
{
    std::ifstream input("../data/"+filename);
    json j;
    input >> j;
    input.close();

    json output_json;
    json children_json;

    for(auto &co:j["CityObjects"].items()) {
        if(co.key() == building_num)
        {
            // information of parents building
            output_json["CityObjects"][building_num] = co.value();
            // cout<<output_json<<endl;
            children_json = co.value()["children"];
        }
    }

    for(auto &children:children_json){
        for(auto &co:j["CityObjects"].items()) {
            if(co.key() == children)
            {
                string name = children;
                json new_j = co.value();
                output_json["CityObjects"][name] = new_j;
            }
        }
    }

    vector<int> old_ver_index;
    //vector<vector<int>> new_vertices;
    int i = 0;

    // for multisurface condition
    for(auto & co:output_json["CityObjects"].items()) {
        for (auto& g : co.value()["geometry"]) {
            if(g["type"] == "MultiSurface") {
                for(auto & boundaries: g["boundaries"]) {
                    for(auto &shells: boundaries) {
                        for(int l = 0; l<shells.size();l++) {

                            if(!old_ver_index.empty()){
                                bool exist = false;
                                for(int k=0;k<old_ver_index.size();k++) {
                                    if(shells[l] == old_ver_index[k]) {

                                        shells[l] = k;
                                        exist = true;
                                        break;
                                    }
                                }
                                if(!exist)
                                {
                                    old_ver_index.emplace_back(shells[l]);
                                    shells[l]=i;
                                    i++;
                                }
                            }
                            else if(old_ver_index.empty()) {
                                old_ver_index.emplace_back(shells[l]);
                                shells[l] = i;
                                i++;
                            }
                        }
                    }
                }
            }
            else if(g["type"] == "Solid") {
                for(auto & boundaries: g["boundaries"]) {
                    for(auto &boundary: boundaries) {
                        for(auto &shells: boundary) {
                            for(int l = 0; l<shells.size();l++) {

                                if(!old_ver_index.empty()){
                                    bool exist = false;
                                    for(int k=0;k<old_ver_index.size();k++) {
                                        if(shells[l] == old_ver_index[k]) {

                                            shells[l] = k;
                                            exist = true;
                                            break;
                                        }
                                    }
                                    if(!exist)
                                    {
                                        old_ver_index.emplace_back(shells[l]);
                                        shells[l]=i;
                                        i++;
                                    }
                                }
                                else if(old_ver_index.empty()) {
                                    cout<<shells<<endl;
                                    old_ver_index.emplace_back(shells[l]);

                                    shells[l] = i;
                                    //new_ver_index.emplace_back(l);
                                    i++;
                                }
                            }
                        }

                    }
                }
            }
        }
    }

    json new_vertices;
    for(auto old_v:old_ver_index) {
        new_vertices.emplace_back(j["vertices"][old_v]);
        //<<new_vertices<<endl;
        //new_vertices.emplace_back(vertices_list[old_v]);
    }

    output_json["metadata"] = j["metadata"];
    output_json["transform"] = j["transform"];
    output_json["type"] = j["type"];
    output_json["version"] = j["version"];
    output_json["vertices"] = new_vertices;

    ofstream o("../data/"+ building_num +".json");
    o << output_json.dump(2) << endl;
    o.close();
    cout<<"finished!"<<endl;
}

json lod_filter(json &j,float lod) {
    
    for (auto &co: j["CityObjects"].items()) {
        json geometry_json;
        for (auto &g: co.value()["geometry"]) {
            if (g["type"] == "Solid") {
                if (g["lod"] == lod)
                    geometry_json.emplace_back(g);
            }
        }
        co.value()["geometry"] = geometry_json;
    }
    return j;
}
// to select surfaces
// in this case, the surface is
static void read_surface_list(json &j, int surface_id)
{
    for(auto &co:j["CityObjects"].items()) {
        for(auto &g:co.value()["geometry"]) {
            if(g["type"] == "Solid") {
                json surfaces_list = g["semantics"]["values"];
                for(auto &surface:surfaces_list){
                }
            }
        }
    }
}


void split_surface(json &j ){

    for (auto& co : j["CityObjects"].items()) {
        for (auto &g: co.value()["geometry"]) {

            if (g["type"] == "Solid") {
                int roof_index = 0;
                int surface_index = 0;
                for (int i = 0; i < g["semantics"]["surfaces"].size(); ++i) {
                    if (g["semantics"]["surfaces"][i]["type"] == "WallSurface")
                    {
                        roof_index = i;
                        surface_index++;
                    }
                    else surface_index++;
                }
                for (auto &values: g["semantics"]["values"]) {

                    for (int k = 0; k < values.size(); ++k) {
                        if (values[k] == roof_index) {
                            for (auto &shell: g["boundaries"]) {
                                g["semantics"]["surfaces"][surface_index]["type"] = "WallSurface";
                                g["semantics"]["surfaces"][surface_index]["area"] = 0.0;
                                g["semantics"]["surfaces"][surface_index]["BuildingPart_id"] = co.key();
                                values[k] = surface_index;
                                surface_index++;
                            }
                        }
                    }
                }
            }
        }
    }
}

// CityJSON files have their vertices compressed: https://www.cityjson.org/specs/1.1.1/#transform-object
// this function visits all the surfaces and print the (x,y,z) coordinates of each vertex encountered
void list_all_vertices(json& j) {
    for (auto& co : j["CityObjects"].items()) {
        std::cout << "= CityObject: " << co.key() << std::endl;
        for (auto& g : co.value()["geometry"]) {
            if (g["type"] == "Solid") {
                for (auto& shell : g["boundaries"]) {
                    for (auto& surface : shell) {
                        for (auto& ring : surface) {
                            std::cout << "---" << std::endl;
                            for (auto& v : ring) {
                                std::vector<int> vi = j["vertices"][v.get<int>()];
                                double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
                                double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
                                double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                                std::cout << std::setprecision(2) << std::fixed << v << " (" << x << ", " << y << ", " << z << ")" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}


// Visit every 'RoofSurface' in the CityJSON model and output its geometry (the arrays of indices)
// Useful to learn to visit the geometry boundaries and at the same time check their semantics.
void visit_roofsurfaces(json &j) {
    for (auto& co : j["CityObjects"].items()) {
        for (auto& g : co.value()["geometry"]) {
            if (g["type"] == "Solid") {
                for (int i = 0; i < g["boundaries"].size(); i++) {
                    for (int j = 0; j < g["boundaries"][i].size(); j++) {
                        int sem_index = g["semantics"]["values"][i][j];
                        if (g["semantics"]["surfaces"][sem_index]["type"].get<std::string>().compare("RoofSurface") == 0) {
                            std::cout << "RoofSurface: " << g["boundaries"][i][j] << std::endl;
                        }
                    }
                }
            }
        }
    }
}

using namespace std;
typedef std::unordered_map<int,Vertex*> vertexmap;
typedef std::map<vector<int>,vector<HalfEdge*>> edgemap;

// 1.
void importOBJ(DCEL & D, const char *file_in) {
    vertexmap id2vertex;
    edgemap ver2halfedge;
    string line;
    char c;
    unsigned int i, j, k;
    float x, y, z;
    string si, sj, sk;
    ifstream in;
    in.open(file_in);
    int n = 0;
    while (getline(in, line))                           // read whole line
    {
        if (line.find_first_of("vf") == std::string::npos) continue;     // skip pointless lines
        istringstream ss(line);                            // put line into a stream for input
        ss >> c;                                             // get first character
        switch (c) {
            case 'v':                                         // VERTICES
            {
                ++n;
                ss >> x >> y >> z;
                id2vertex[n] = D.createVertex(x, y, z);
                D.vermap()[id2vertex[n]] = n - 1;
                break;
            }
            case 'f':                                         // FACES
            {
                Face *f0 = D.createFace();
                ss >> si >> sj >> sk;
                i = std::stoi(si);
                j = std::stoi(sj);
                k = std::stoi(sk);
                HalfEdge *e0 = D.createHalfEdge();
                HalfEdge *e1 = D.createHalfEdge();
                HalfEdge *e2 = D.createHalfEdge();
                e0->origin = id2vertex[i];
                e0->destination = id2vertex[j];
                vector<int> l0;
                l0.push_back(i);
                l0.push_back(j);
                sort(l0.begin(),l0.end());
                ver2halfedge[l0].push_back(e0);
                e0->prev = e2;
                e0->next = e1;
                f0->exteriorEdge = e0;
                e0->incidentFace = f0;

                e1->origin = id2vertex[j];
                e1->destination = id2vertex[k];
                vector<int> l1;
                l1.push_back(j);
                l1.push_back(k);
                sort(l1.begin(),l1.end());
                ver2halfedge[l1].push_back(e1);
                e1->prev = e0;
                e1->next = e2;
                e1->incidentFace = f0;

                e2->origin = id2vertex[k];
                e2->destination = id2vertex[i];
                vector<int> l2;
                l2.push_back(k);
                l2.push_back(i);
                sort(l2.begin(),l2.end());
                ver2halfedge[l2].push_back(e2);
                e2->prev = e1;
                e2->next = e0;
                e2->incidentFace = f0;
                break;
            }
        }
    }
    in.close();

    edgemap::iterator iter;
    iter = ver2halfedge.begin();
    for (iter = ver2halfedge.begin(); iter != ver2halfedge.end(); iter++) {
        auto ver_list = iter->first;
        vector<HalfEdge*> edge_list = iter->second;

        if (edge_list.size() != 2) throw std::invalid_argument("Non-manifold!!!");
        HalfEdge* edge1 = edge_list[0];
        HalfEdge* edge2 = edge_list[1];
        edge1->twin = edge2;
        edge2->twin = edge1;

//        if(edge_list.size() !=2)
//        {
//            cout<<"error"<<endl;
//            cout<<ver_list[0]<<", "<<ver_list[1]<<endl;
//        }
//        else
//        {
//            HalfEdge* edge1 = edge_list[0];
//            HalfEdge* edge2 = edge_list[1];
//            edge1->twin = edge2;
//            edge2->twin = edge1;
//        }
    }
    //printDCEL(D);
}

// 2.
void groupTriangles(DCEL & D) {
    int n = 0;
    int i;
    list<Face*> facelist_notest;
    facelist_notest = (const list<Face*, allocator<Face*>> &) D.faces();
    while (facelist_notest.size() > 0) {
        n++;
        std::vector<Face*> check_list;
        auto iter = facelist_notest.begin();
        Face* test_face = *iter;
        check_list.push_back(test_face);
        facelist_notest.remove(test_face);
        test_face->meshnum = n;
        D.meshes()[n].push_back(test_face);
        for (i = 0; i < check_list.size(); i++) {
            HalfEdge* e = check_list[i]->exteriorEdge;
            const HalfEdge* e_start = e;
            do {
                if (((e->twin)->incidentFace)->meshnum == n)
                    e = e->next;
                else {
                    check_list.push_back((e->twin)->incidentFace);
                    facelist_notest.remove((e->twin)->incidentFace);
                    ((e->twin)->incidentFace)->meshnum = n;
                    D.meshes()[n].push_back((e->twin)->incidentFace);
                    e = e->next;
                }
            } while (e_start != e);
        }
    }
}

// 3.
float signed_volume(const Point& a, const Point& b, const Point& c, const Point& d) {
    Point a_d = a - d;
    Point b_d = b - d;
    Point c_d = c - d;
    return a_d.dot(b_d.cross(c_d));
}

Face* flipOrient(Face* face) {
    Vertex* v0 = face->exteriorEdge->origin;
    Vertex* v1 = face->exteriorEdge->destination;
    Vertex* v2 = face->exteriorEdge->next->destination;
    HalfEdge* e0 = face->exteriorEdge;
    HalfEdge* e1 = face->exteriorEdge->next;
    HalfEdge* e2 = face->exteriorEdge->prev;
    e0->origin = v1;
    e0->destination = v0;
    e0->next = e2;
    e0->prev = e1;
    e1->origin = v2;
    e1->destination = v1;
    e1->next = e0;
    e1->prev = e2;
    e2->origin = v0;
    e2->destination = v2;
    e2->next = e1;
    e2->prev = e0;
    return face;
}

void orientMeshes(DCEL & D) {
    for (auto iter = D.meshes().begin(); iter != D.meshes().end(); ++iter) {
        list<Face*> value = iter->second;
        Vertex* v_highest = ((*value.begin())->exteriorEdge)->destination;
        for (auto iter_f = value.begin(); iter_f != value.end(); ++iter_f) {
            HalfEdge* e = (*iter_f)->exteriorEdge;
            const HalfEdge* e_start = e;
            Vertex* v_high = e->destination;
            auto z_max = v_high ->z;
            do {
                e = e->next;
                if ((e->destination)->z > z_max)
                    v_high = e->destination;
                if (v_high->z > v_highest->z)
                    v_highest = v_high;
            } while (e_start != e);
        }
        const HalfEdge* e_slope = (*value.begin())->exteriorEdge;
        float m = 1;
        for (auto iter_f = value.begin(); iter_f != value.end(); ++iter_f) {
            HalfEdge* e = (*iter_f)->exteriorEdge;
            const HalfEdge* e_start = e;
            do {
                if ((e->destination) == v_highest) {
                    if (((v_highest->z) - ((e->origin)->z)) / sqrt(((v_highest->x) - ((e->origin)->x)) * ((v_highest->x) - ((e->origin)->x)) + ((v_highest->y) - ((e->origin)->y)) * ((v_highest->y) - ((e->origin)->y))) < m){
                        m = ((v_highest->z) - ((e->origin)->z)) / sqrt(((v_highest->x) - ((e->origin)->x)) * ((v_highest->x) - ((e->origin)->x)) + ((v_highest->y) - ((e->origin)->y)) * ((v_highest->y) - ((e->origin)->y)));
                        e_slope = e;
                    }
                }
                if ((e->origin) == v_highest) {
                    if (((v_highest->z) - ((e->destination)->z)) / sqrt(((v_highest->x) - ((e->destination)->x)) * ((v_highest->x) - ((e->destination)->x)) + ((v_highest->y) - ((e->destination)->y)) * ((v_highest->y) - ((e->destination)->y))) < m) {
                        m = ((v_highest->z) - ((e->destination)->z)) / sqrt(((v_highest->x) - ((e->destination)->x)) * ((v_highest->x) - ((e->destination)->x)) + ((v_highest->y) - ((e->destination)->y)) * ((v_highest->y) - ((e->destination)->y)));
                        e_slope = e;
                    }
                }
                e = e->next;
            } while (e_start != e);
        }

        Face* f_intersect = e_slope->incidentFace;

        if (((((e_slope->twin)->next)->destination)->z) > (((e_slope->next)->destination)->z))
            f_intersect = (e_slope->twin)->incidentFace;

        Vertex* v0 = f_intersect->exteriorEdge->origin;
        Vertex* v1 = f_intersect->exteriorEdge->destination;
        Vertex* v2 = f_intersect->exteriorEdge->next -> destination;
        Point p0(v0->x, v0->y, v0->z);
        Point p1(v1->x, v1->y, v1->z);
        Point p2(v2->x, v2->y, v2->z);
        Point p3(v_highest->x, v_highest->y, (v_highest->z) + 10.0f);
        if (signed_volume(p0, p1, p2, p3) > 0) {
            flipOrient(f_intersect);
        }
        if (signed_volume(p0, p1, p2, p3) == 0) {
            f_intersect = e_slope->incidentFace;
            Vertex* v0 = f_intersect->exteriorEdge->origin;
            Vertex* v1 = f_intersect->exteriorEdge->destination;
            Vertex* v2 = f_intersect->exteriorEdge->next->destination;
            Point p0(v0->x, v0->y, v0->z);
            Point p1(v1->x, v1->y, v1->z);
            Point p2(v2->x, v2->y, v2->z);
            Point p3(v_highest->x, v_highest->y, (v_highest->z) + 10.0f);
            if (signed_volume(p0, p1, p2, p3) > 0) {
                flipOrient(f_intersect);
            }
        }
        f_intersect->flip = 1;
        std::vector<Face*> check_list;
        check_list.push_back(f_intersect);
        int i;
        for (i = 0; i < check_list.size(); i++) {
            HalfEdge* e = check_list[i]->exteriorEdge;
            const HalfEdge* e_start = e;
            do {
                if (((e->twin)->origin) == e->destination && e->twin->incidentFace->flip == 1) {
                    e = e->next;
                }
                else if (((e->twin)->origin) == e->destination && e->twin->incidentFace->flip == 0) {
                    e->twin->incidentFace->flip = 1;
                    check_list.push_back(e->twin->incidentFace);
                    e = e->next;
                }
                else if (((e->twin)->origin) != e->destination) {
                    e->twin->incidentFace->flip = 1;
                    flipOrient(e->twin->incidentFace);
                    check_list.push_back(e->twin->incidentFace);
                    e = e->next;
                }

            } while (e_start != e);
        }
    }
}

// 4.
Point get_Normal(Face* & f) {
    const HalfEdge* e = f->exteriorEdge;
    const Vertex * v0 = e->origin;
    const Vertex * v1 = e->destination;
    const Vertex * v2 = (e->next)->destination;
    float a = ((v1->y - v0->y) * (v2->z - v0->z) - (v1->z - v0->z) * (v2->y - v0->y));
    float b = ((v1->z - v0->z) * (v2->x - v0->x) - (v1->x - v0->x) * (v2->z - v0->z));
    float c = ((v1->x - v0->x) * (v2->y - v0->y) - (v1->y - v0->y) * (v2->x - v0->x));
    //normalized
    float length = sqrt(a*a + b*b + c*c);
    a = a / length;
    b = b / length;
    c = c / length;
    Point norm_vec{ a,b,c };
    return norm_vec;
}

bool isBoundary (HalfEdge* & edge, const Point &this_normal, list<HalfEdge*> & boundaryList){
    //check if an edge is a boundary edge for a co-planar face
    Face *face1 = edge->twin->incidentFace;
    Point normal1 = get_Normal(face1);
    bool isboundary = false;
    Point cross_vec = this_normal.cross(normal1);
    float cross_length = this_normal.cross_length(normal1);
//    if(!(this_normal == normal1 || ((fabs(cross_vec.x) <= 1e-5) && (fabs(cross_vec.y) <= 1e-5) && (fabs(cross_vec.z) <= 1e-5)))){
    if(!(this_normal == normal1 || cross_length <= 1e-1 || ((fabs(cross_vec.x) <= 1e-5) && (fabs(cross_vec.y) <= 1e-5) && (fabs(cross_vec.z) <= 1e-5)))){
        boundaryList.push_back(edge);
        isboundary = true;
    }
    return isboundary;
}

void  coPlanarTest (list<Face*> & mesh_value_list, std::list<Face*>& test_list, Face* & test_face,list<HalfEdge*> & boundaryList) {
    std::queue<HalfEdge*> check_queue;
    HalfEdge *e1 = test_face->exteriorEdge;
    HalfEdge *e2 = e1->next;
    HalfEdge *e3 = e1->prev;
    Point this_normal = get_Normal(test_face);
    check_queue.push(e1);
    check_queue.push(e2);
    check_queue.push(e3);
    int n = 0;
    while(n < 3){
        ++n;
        HalfEdge *e0 = check_queue.front();
        check_queue.pop();
        if(isBoundary(e0, this_normal, boundaryList)){
            e0->incidentFace = test_face;
        }
        else{
            test_list.remove(e0->twin->incidentFace);
            mesh_value_list.remove(e0->twin->incidentFace);
            e0->eliminate();
            e0->twin->eliminate();
            e0->next->prev = e0->twin->prev;
            e0->prev->next = e0->twin->next;
            e0->twin->next->prev = e0->prev;
            e0->twin->prev->next = e0->next;
            (e0->twin)->incidentFace->eliminate();
            (e0->twin)->next->incidentFace = test_face; // it is ok to delete?
            (e0->twin)->prev->incidentFace = test_face;
            check_queue.push(e0->twin->next);
            check_queue.push(e0->twin->prev);
        }
    }

    while (!check_queue.empty()) {
//        std::cout << "test check_queue.size()= "<< check_queue.size() << std::endl;
        HalfEdge* e = check_queue.front();
        Face* tf = e->twin->incidentFace;
        check_queue.pop();
        if (isBoundary(e, this_normal, boundaryList)){ e->incidentFace = test_face; }
        else if (tf != test_face){
            test_list.remove(tf);
            mesh_value_list.remove(tf);
            e->eliminate();
            e->twin->eliminate();
            tf->eliminate();
            e->next->prev = e->twin->prev;
            e->prev->next = e->twin->next;
            e->twin->next->prev = e->prev;
            e->twin->prev->next = e->next;
            (e->twin)->next->incidentFace = test_face;
            (e->twin)->prev->incidentFace = test_face;
            e->next->incidentFace = test_face;
            e->prev->incidentFace = test_face;
            check_queue.push(e->twin->next);
            check_queue.push(e->twin->prev);
        }
        else if ( tf == test_face){
            e->eliminate();
            e->twin->eliminate();
            e->next->prev = e->twin->prev;
            e->prev->next = e->twin->next;
            e->twin->next->prev = e->prev;
            e->twin->prev->next = e->next;
        }
    }
//        do {
//            Face* & tf = e->twin->incidentFace;
//            if(isBoundary(e, this_normal, boundaryList)){
//                e->incidentFace = test_face;
//                e = e->next;}
//            else if (e->isEliminated()){e = e->next;}
//            else if (tf != test_face){
//                check_list.push_back(tf);
//                test_list.remove(tf);
//                mesh_value_list.remove(tf);
//                e->next->prev = e->twin->prev;
//                e->prev->next = e->twin->next;
//                e->twin->next->prev = e->prev;
//                e->twin->prev->next = e->next;
//                e->eliminate();
//                e->twin->eliminate();
//                (e->twin)->incidentFace->eliminate();
//                (e->twin)->incidentFace = test_face;
//                (e->twin)->next->incidentFace = test_face; // it is ok to delete
//                (e->twin)->prev->incidentFace = test_face;
//                e->incidentFace = test_face;
//                e->next->incidentFace = test_face;
//                e->prev->incidentFace = test_face;
//                e = e->next;
//            }
//            else if (tf == test_face){
//                e->eliminate();
//                e->twin->eliminate();
//                e->next->prev = e->twin->prev;
//                e->prev->next = e->twin->next;
//                e->twin->next->prev = e->prev;
//                e->twin->prev->next = e->next;
//                e = e->next;
//            }
//            else{e = e-> next;}
//        } while ( e_start!=e) ;
}

void computeBbox(const std::vector<HalfEdge*> & edgelist,
                 double& max_x, double& max_y, double& max_z, double& min_x, double& min_y, double& min_z){
    min_x = max_x = edgelist[0]->origin->x;
    min_y = max_y = edgelist[0]->origin->y;
    min_z = max_z = edgelist[0]->origin->z;
    for (auto &each: edgelist){
        Vertex* v0 = each->origin;
        Vertex* v1 = each->destination;
        if(v0->x > max_x) max_x = v0->x;
        else if (v0->x < min_x) min_x = v0->x;
        if(v0->y > max_y) max_y = v0->y;
        else if (v0->y < min_y) min_y = v0->y;
        if(v0->z > max_z) max_z = v0->z;
        else if (v0->z < min_z) min_z = v0->z;
        if(v1->x > max_x) max_x = v1->x;
        else if (v1->x < min_x) min_x = v1->x;
        if(v1->y > max_y) max_y = v1->y;
        else if (v1->y < min_y) min_y = v1->y;
        if(v1->z > max_z) max_z = v1->z;
        else if (v1->z < min_z) min_z = v1->z;
    }
//    std::cout << "test xyz: " << min_x << ", " << min_y << ", " << min_z
//    << ", " << max_x << ", " << max_y << ", " << max_z << std::endl;
}

void mergeCoPlanarFaces(DCEL & D) {
    for (auto iter = D.meshes().begin(); iter != D.meshes().end(); ++iter) {
        list<Face*> & value = iter->second;
        list<Face*> test_list = iter->second;
        while(!test_list.empty()){
            list<HalfEdge*> boundaryList;
            Face* test_face = *(test_list.begin());
            test_list.remove(test_face);

            coPlanarTest(value, test_list,test_face, boundaryList);

            std::vector<std::vector<HalfEdge*>> boundary_hole_List;//a list include exterior boundary and holes
            if (boundaryList.size() > 3){
                int n = 0;
                while(!boundaryList.empty()){
                    auto iter = boundaryList.begin();
                    std::vector<HalfEdge*> this_list;
                    this_list.push_back(*iter);
                    boundary_hole_List.push_back(this_list);
                    HalfEdge* e = *iter;
                    HalfEdge* e_start = e;
                    boundaryList.remove(*iter);
                    do{
                        boundary_hole_List[n].push_back(e->next);
                        boundaryList.remove(e->next);
                        e = e->next;
                    }
                    while(e_start != e);
                    ++n;
                }
            }
//            else{throw std::invalid_argument("A face with only a single triangle");}

            //check how many holes a plane has
//            std::cout << "test boundary_hole_List.size()= " << boundary_hole_List.size() << std::endl;
            if(boundary_hole_List.size() > 1){
                double max_x, max_y, max_z, min_x, min_y, min_z;
                computeBbox(boundary_hole_List[0],max_x, max_y, max_z, min_x, min_y, min_z);
                test_face->holes = {};
                for (int n = 1; n < boundary_hole_List.size(); ++n){
                    double temp_max_x, temp_max_y, temp_max_z, temp_min_x, temp_min_y, temp_min_z;
                    computeBbox(boundary_hole_List[n],temp_max_x, temp_max_y, temp_max_z, temp_min_x, temp_min_y, temp_min_z);
                    if (temp_max_x >= max_x && temp_max_y >= max_y && temp_max_z >= max_z &&
                        temp_min_x <= min_x && temp_min_y <= min_y && temp_min_z <= min_z){
                        test_face->exteriorEdge = boundary_hole_List[n][0];
                        test_face->holes.push_back(boundary_hole_List[0][0]);
                    }
                    else if (temp_max_x <= max_x && temp_max_y <= max_y && temp_max_z <= max_z &&
                             temp_min_x >= min_x && temp_min_y >= min_y && temp_min_z >= min_z){
                        test_face->exteriorEdge = boundary_hole_List[0][0];
                        test_face->holes.push_back(boundary_hole_List[n][0]);
                    }
                    else{
                        test_face->holes.push_back(boundary_hole_List[0][0]);
                        test_face->holes.push_back(boundary_hole_List[n][0]);
                    }
                }
            }
            else if (boundary_hole_List.size() == 1){
                test_face->exteriorEdge = boundary_hole_List[0][0];
                for (auto each: boundary_hole_List[0]){
                    each->incidentFace = test_face;
                }
            }
        }
    }
    D.cleanup();
//    printDCEL(D);
//    for (const auto &each: D.halfEdges()){
//        if (each->hasDanglingLink()){
//            std::cout << std::boolalpha;
//            std::cout << "test danglinglink: " << each->isEliminated() << std::endl;
//            std::cout << "test danglinglink: " << each->incidentFace->isEliminated() << ", " << each->origin->isEliminated() << ", " << each->destination->isEliminated()
//                      << ", " << each->next->isEliminated() << ", " << each->prev->isEliminated() << ", " << each->twin ->isEliminated()<< std::endl;
//        }
//    }
}

// 5.
void exportCityJSON(DCEL& D, const char* file_out) {

    int num_hole = 0;

    std::ofstream outfile(file_out, std::ios::out);
    int i;
    outfile << "{" << std::endl
            << "\"type\":\"CityJSON\"," << std::endl
            << "\"version\":\"1.0\"," << std::endl
            << "\"CityObjects\":{" << std::endl;

    for (auto iter = D.meshes().begin(); iter != D.meshes().end(); ++iter) {
        int key = iter->first;
        list<Face*> value = iter->second;

        outfile << "\"id-" << key << "\": {" << std::endl
                << "\"type\": \"BuildingPart\"," << std::endl
                << "\"geometry\": [{" << std::endl
                << "\"type\": \"MultiSurface\"," << std::endl
                << "\"lod\": 2.0," << std::endl
                << "\"boundaries\": [" << std::endl;
        auto iter_v = value.begin();
        HalfEdge* e = (*iter_v) ->exteriorEdge;
        const HalfEdge* e_start = e;

        outfile << "[[" << D.vermap()[e->origin];

        do {
            e = e->next;
            outfile << ", " << D.vermap()[e->origin];
        } while (e_start != e->next);
        outfile << "]";

        if (((*iter_v)->holes).empty())
            outfile << "]";
        else {
            for (auto iter_hole = (*iter_v)->holes.begin(); iter_hole != (*iter_v)->holes.end(); ++iter_hole) {

                ++num_hole;

                outfile << ", [";
                HalfEdge* e = *iter_hole;
                const HalfEdge* e_start = e;
                outfile <<  D.vermap()[e->origin];
                do {
                    e = e->next;
                    outfile << ", " << D.vermap()[e->origin];
                } while (e_start != e->next);
                outfile << "]";
            }
            outfile << "]";
        }
        ++iter_v;
        while (iter_v != value.end()) {
            HalfEdge* e = (*iter_v)->exteriorEdge;
            const HalfEdge* e_start = e;
            outfile << ", [[" << D.vermap()[e->origin];
            do {
                e = e->next;
                outfile << ", " << D.vermap()[e->origin];
            } while (e_start != e->next);
            outfile << "]";

            if (((*iter_v)->holes).empty())
                outfile << "]";
            else {
                for (auto iter_hole = (*iter_v)->holes.begin(); iter_hole != (*iter_v)->holes.end(); ++iter_hole) {

                    ++num_hole;

                    outfile << ", [";
                    HalfEdge* e = *iter_hole;
                    const HalfEdge* e_start = e;
                    outfile << D.vermap()[e->origin];
                    do {
                        e = e->next;
                        outfile << ", " << D.vermap()[e->origin];
                    } while (e_start != e->next);
                    outfile << "]";
                }
                outfile << "]";
            }

            ++iter_v;
        }

        auto tmp = iter;
        advance(tmp,1);

        if (tmp != D.meshes().end())
            outfile << "]" << std::endl << "}]" << std::endl << "}," << std::endl;
        else
            outfile << "]" << std::endl << "}]" << std::endl << "}" << std::endl;

//        if(iter == D.meshes().end())
//            outfile << "]" << std::endl << "}]" << std::endl << "}" << std::endl;
//        else
//            outfile << "]" << std::endl << "}]" << std::endl << "}," << std::endl;
    }
    outfile << "}," << std::endl;
    outfile << "\"vertices\": [" << std::endl;
    auto iter = D.vertices().begin();
    outfile << "[" << (*iter)->x << ", " << (*iter)->y << ", " << (*iter)->z << "]";
    ++iter;
    while(iter != D.vertices().end()) {
        outfile << ", [" << (*iter)->x << ", " << (*iter)->y << ", " << (*iter)->z << "]";
        ++iter;
    }
    outfile << std::endl;
    outfile << "]" << std::endl;
    outfile << "}" << std::endl;

    outfile.close();
//std::cout << "total holes: " << num_hole << std::endl;
}



int main(int argc, const char * argv[]) {

    string filename = "3dbag_v210908_fd2cee53_5910.json";
    string building_num = "NL.IMBAG.Pand.0503100000032914";

    std::ifstream input("../data/"+filename);
    json j;
    input >> j;
    input.close();
    
    lod_filter(j,2.2);
    select_single_building(filename,building_num);

    clock_t start, end;
    start = clock();
    const char *file_in = "..//data//NL.IMBAG.Pand.0503100000018507_lod22_tri.obj";
    const char *file_out = "..//data//result.json";

    // create an empty DCEL
    DCEL D;
    // 1. read the triangle soup from the OBJ input file and convert it to the DCEL,
    importOBJ(D, file_in);
    // 2. group the triangles into meshes,
    groupTriangles(D);
    // 3. determine the correct orientation for each mesh and ensure all its triangles
    //    are consistent with this correct orientation (ie. all the triangle normals
    //    are pointing outwards).
    orientMeshes(D);

    // 4. merge adjacent triangles that are co-planar into larger polygonal faces.
    mergeCoPlanarFaces(D);
    // 5. write the meshes with their faces to a valid CityJSON output file.
    exportCityJSON(D, file_out);

    end = clock();
    std::cout << "Time " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
    
    
    return 0;
}
