#include <iostream>
#include <fstream>
#include <string>
#include "include/json.hpp"
#include <cmath>


using namespace std;
using json = nlohmann::json;

class Point{
public:
    double x;
    double y;
    double z;

    Point(double x1,double y1,double z1): x(x1),y(y1),z(z1)
    {}
};

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




int main(int argc, const char * argv[]) {

    string filename = "3dbag_v210908_fd2cee53_5910.json";
    string building_num = "NL.IMBAG.Pand.0503100000032914";

    std::ifstream input("../data/"+filename);
    json j;
    input >> j;
    input.close();
    
    lod_filter(j,2.2);
    select_single_building(filename,building_num);
    
    
    return 0;
}
