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
#include "include/json.hpp"
#include <CGAL/Vector_3.h>

#include "DCEL.hpp"


using namespace std;
using json = nlohmann::json;

#define M_PI  3.141592653589793238462643383279502884L;

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

float distance(Point &pt1, Point &pt2)
{
    float dis = sqrt(pow((pt1.x-pt2.x),2)+pow((pt1.y-pt2.y),2)+pow((pt1.z-pt2.z),2));
    return dis;
}

//
int nearest_surface(json &j, Point &camera) {
    // find out windows

    vector<Point> centre;
    for (auto& co : j["CityObjects"].items()) {
        std::cout << "= CityObject: " << co.key() << std::endl;
        for (auto& g : co.value()["geometry"]) {
            if (g["type"] == "MultiSurface") {
                for (auto& shell : g["boundaries"]) {
                    Point pt(0.0,0.0,0.0);
                    int count=0;
                    for(auto& surface: shell)
                    {
                        for (auto& ring : surface) {
                            for (auto& v : ring) {
                                std::vector<int> vi = j["vertices"][v.get<int>()];
                                pt.x += vi[0];
                                pt.y += vi[1];
                                pt.z += vi[2];
                                count = count+1;
                            }
                        }
                    }
                    centre.emplace_back(pt.x/count,pt.y/count,pt.z/count);
                }
            }
        }
    }
    int count=0;
    float nearest_distance= distance(camera,centre[0]);
    for(int i=0;i<centre.size();i++) {
        float current_dis = distance(camera,centre[i]);
        cout<<current_dis<<endl;
        if(nearest_distance >= current_dis)
        {
            nearest_distance = current_dis;
            count = i;
        }
    }
    cout<<count<<endl;
    return count;
}

// function to process window corners
//output: a vector<vector<double*>>, the order of corner point is:
//top-left, bottom-left, bottom-right, top-right, in CCW.
vector<vector<Point>> window_processing(string &filename){
    string line;
    ifstream in;
    in.open(filename);
    vector<Point> vertices;
    vector<vector<Point>> window_list;
    if(in.is_open())
    {
        std::string line;
        while (getline(in, line)) {
            std::istringstream iss(line);
            std::string word;
            iss >> word;
            vector<float> coordinates;
            if (word == "window" || word == "door") {
                if(vertices.empty()) continue;
                else{
                    window_list.emplace_back(vertices);
                    vertices.clear();
                }
            }
            else {
                coordinates.push_back(std::stof(word));
                while (iss >> word)
                    coordinates.push_back(std::stof(word));
                if (coordinates.size() == 3) vertices.emplace_back(coordinates[0], coordinates[1], coordinates[2]);
                else vertices.emplace_back();
            }
        }
    }
    if(!vertices.empty())
        window_list.emplace_back(vertices);
    return window_list;
}

json add_multi_window(string &building_file, string &window_file, Point &camera)
{
    ifstream input("../data/"+ building_file);
    json j;
    input >> j;
    input.close();

//    int surface_no = nearest_surface(j,camera);
    vector<vector<Point>> objects_corner = window_processing(window_file);

    // add to vertices list

    for(auto &window_corner:objects_corner)
    {
        json window_index;
        int size = j["vertices"].size();
        for(int i = 0;i < window_corner.size(); i++)
        {
            json vertex = {window_corner[i].x,window_corner[i].y,window_corner[i].z};
            cout<<vertex<<endl;
            //j["vertices"].emplace_back(vertex);
            j["vertices"][size+i]={window_corner[i].x,window_corner[i].y,window_corner[i].z};
            cout<<j["vertices"][size+i]<<endl;
            window_index[i]=size+i;
        }

        for (auto& co : j["CityObjects"].items()) {
            std::cout << "= CityObject: " << co.key() << std::endl;
            for (auto& g : co.value()["geometry"]) {
                if (g["type"] == "MultiSurface") {
                    json tmp = {window_index};
                    g["boundaries"].emplace_back(tmp);
                }
            }
        }
    }

    return j;
}


void calculate_centroid(string &building_file)
{

    ifstream input("../data/"+building_file);
    json j;
    input >> j;
    input.close();

    vector<Point> centre;
    string building_name;

    for (auto& co : j["CityObjects"].items()) {
        //std::cout << "= CityObject: " << co.key() << std::endl;
        for (auto &g: co.value()["geometry"]) {
            if (g["lod"] == 1.2)
            {
                Point bbox_max(0,0,0);
                Point bbox_min(999999,999999,0);
                int count=0;
                for (auto &shell:g["boundaries"]) {
                    for(auto &item:shell){
                        for(auto & surfaces:shell)
                        {
                            for(auto &surface:surfaces) {
                                for (auto& v : surface) {
                                    Point pt(0.0,0.0,0.0);
                                    std::cout<<v<<std::endl;
                                    vector<int> vi = j["vertices"][v.get<int>()];
                                    pt.x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
                                    pt.y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
                                    pt.z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                                    count = count+1;
                                    if(bbox_max.x<pt.x) bbox_max.x = pt.x;
                                    if(bbox_max.y<pt.y) bbox_max.y = pt.y;
                                    if(bbox_min.x>pt.x) bbox_min.x = pt.x;
                                    if(bbox_min.y>pt.y) bbox_min.y = pt.y;
                                }
                            }
                        }
                    }
                }
                Point center ((bbox_max.x+bbox_min.x)/2,(bbox_max.y+bbox_min.y)/2,2);
                cout<<"the center coordinates is "<<center.x<<", "<<center.y<<", "<<endl;

                std::ofstream outfile("centroid.json", std::ios::out);

                }
            }
        }
    }

// for 3d buildings, lod 0 is not always correct.
// we might need the 2d bbox of the building
string buildings_determination(string &building_file, Point &camera, float deg) {

    ifstream input("../data/"+building_file);
    json j;
    input >> j;
    input.close();

    vector<Point> centre;
    string building_name;

    json output_j;
    float nearest_dis=99999;

    for (auto& co : j["CityObjects"].items()) {
        //std::cout << "= CityObject: " << co.key() << std::endl;
        for (auto &g: co.value()["geometry"]) {
            if (g["lod"] == 1.2)
            {
                Point bbox_max(0,0,0);
                Point bbox_min(999999,999999,0);
                int count=0;
                for (auto &shell:g["boundaries"]) {
                    for(auto &item:shell){
                        for(auto & surfaces:shell)
                        {
                            for(auto &surface:surfaces) {
                                for (auto& v : surface) {
                                    Point pt(0.0,0.0,0.0);
                                    //std::cout<<v<<std::endl;
                                    vector<int> vi = j["vertices"][v.get<int>()];
                                    pt.x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
                                    pt.y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
                                    pt.z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                                    count = count+1;
                                    if(bbox_max.x<pt.x) bbox_max.x = pt.x;
                                    if(bbox_max.y<pt.y) bbox_max.y = pt.y;
                                    if(bbox_min.x>pt.x) bbox_min.x = pt.x;
                                    if(bbox_min.y>pt.y) bbox_min.y = pt.y;
                                }
                            }
                        }
                    }
                }

                Point center ((bbox_max.x+bbox_min.x)/2,(bbox_max.y+bbox_min.y)/2,2);
                //cout<<"the center coordinates is "<<center.x<<", "<<center.y<<", "<<endl;
                float delta_x=center.x-camera.x;
                float delta_y=center.y-camera.y;
                float rad = atan(delta_x/delta_y);

                float angle = rad/3.141592653589793238462643383279502884 * 180.0;


                float current_dis= distance(camera,center);

                if(angle <= deg+45 && angle>=deg-45) {

                    if(current_dis<nearest_dis)
                    {
                        nearest_dis = current_dis;
                        building_name=co.key();
                        cout<<"the building is"<<co.key();
                        cout<<"  the angle is: "<<angle;
                        cout<<" the distance is: "<<current_dis<<endl;
                    }

                }
            }
        }

    }
    cout<<(building_name)<<endl;
    return building_name;
}


void surface_determination(json &j,string &building_name) {
    for (auto& co : j["CityObjects"].items()) {
        //std::cout << "= CityObject: " << co.key() << std::endl;
        if(co.key() == building_name)
        {
            
        }
    }
}

class LoD3 {
public:
    json building_j;
    Point camera;
    vector<vector<Point>> windows_corner;
    //json add_multi_window(string &building_file, string &window_file, Point &camera)
    json add_windows(string &building_file, string &window_file, Point &camera);
    //string buildings_determination(string &building_file, Point &camera,float & camera_direction);

};

json LoD3::add_windows(string &building_file, string &window_file, Point &camera)
{
    ifstream input("../data/"+ building_file);
    json j;
    input >> j;
    input.close();

//    int surface_no = nearest_surface(j,camera);
    vector<vector<Point>> objects_corner = window_processing(window_file);

    // add to vertices list

    for(auto &window_corner:objects_corner)
    {
        json window_index;
        int size = j["vertices"].size();
        for(int i = 0;i < window_corner.size(); i++)
        {
            json vertex = {window_corner[i].x,window_corner[i].y,window_corner[i].z};
            cout<<vertex<<endl;
            //j["vertices"].emplace_back(vertex);
            j["vertices"][size+i]={window_corner[i].x,window_corner[i].y,window_corner[i].z};
            cout<<j["vertices"][size+i]<<endl;
            window_index[i]=size+i;
        }

        for (auto& co : j["CityObjects"].items()) {
            std::cout << "= CityObject: " << co.key() << std::endl;
            for (auto& g : co.value()["geometry"]) {
                if (g["type"] == "MultiSurface") {
                    json tmp = {window_index};
                    g["boundaries"].emplace_back(tmp);
                }
            }
        }
    }

    return j;
}



int main(int argc, const char * argv[]) {

    string filename = "5910.json";
    string window_file = "..//data//windows.txt";

    Point camera(	85311.87,446874.23,1);

    // deg could be obtained from xml file
    string building_name = buildings_determination(filename,camera,77.965);

    ifstream input("../data/"+filename);
    json j;
    input >> j;
    input.close();

    json new_j = add_multi_window(filename,window_file,camera);
    ofstream output("..//data//test.json");
    output<<new_j<<endl;
    output.close();


    return 0;
}
